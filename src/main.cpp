#include "io.h"
#include "geom.h"

#include "gitversion.h"

#include <fstream>
#include <iostream>
#include <chrono>
#include <getopt.h>
#include <csignal>

INITIALIZE_EASYLOGGINGPP
unsigned DBG_INDENT_CTR = 0;
std::default_random_engine random_engine;
bool main_loop_interrupted = false;

/*seconds*/

static void
setup_logging(int argc, char* argv[]) {
  START_EASYLOGGINGPP(argc, argv);

  el::Configurations defaultConf;
  defaultConf.setGlobally(el::ConfigurationType::Format, "%datetime{%H%m%s.%g} %levshort %msg");
  el::Loggers::reconfigureAllLoggers(defaultConf);
  el::Loggers::addFlag( el::LoggingFlag::ColoredTerminalOutput );
}

[[noreturn]]
static void
usage(const char *progname, int err) {
  std::ostream &f = err ? std::cerr : std::cout;

  f << "Usage: " << progname << "[options] <INPUT> <OUTPUT>" << std::endl
    << "  Options" << std::endl
    << "    --seed NUM             seed of the RNG" << std::endl
    << "    --full-obj             also print vertex coordinates to .obj file" << std::endl
    << "    --face-obj             outout a face-based obj (instead of a segment based)" << std::endl
    << "    --to-beat     TO_BEAT  return when a better answer is found" << std::endl
    << "    --lower-bound NUM      return immediately when this is reached" << std::endl
    << "    --improve RUNS         Do at least RUNS runs/attempts at solving this and improving it *after* beating TO_BEAT" << std::endl
    << "    --improve-time SECONDS After beating TO_BEAT, work on it for at least this long" << std::endl
    << "    --improve-max RUNS     Try at most RUNS runs/attempts at solving this before beating TO_BEAT" << std::endl
    << "    --max-time NUM         Do not start a new run after NUM seconds (overrides improve-* bounds)" << std::endl
    << "    --log-interval SECONDS Report on state regularly." << std::endl
    << "    --obj-in               Input is an obj file, potentially with already segments/faces to improve" << std::endl
  ;
  exit(err);
}


static void
signalHandler( int signum ) {
   LOG(INFO) << "Interrupt signal (" << signum << ") received.\n";
   main_loop_interrupted = true;
}

int main(int argc, char *argv[]) {
  const char * const short_options = "hS:fb:B:I:i:T:L:M:OF";
  const option long_options[] = {
    { "help"        , no_argument      , 0, 'h'},
    { "seed"        , required_argument, 0, 'S'},
    { "full-obj"    , no_argument      , 0, 'f'},
    { "face-obj"    , no_argument      , 0, 'F'},
    { "to-beat"     , required_argument, 0, 'b'},
    { "lower-bound" , required_argument, 0, 'B'},
    { "improve"     , required_argument, 0, 'I'},
    { "improve-time", required_argument, 0, 'i'},
    { "improve-max" , required_argument, 0, 'M'},
    { "max-time"    , required_argument, 0, 'T'},
    { "log-interval", required_argument, 0, 'L'},
    { "obj-in",       no_argument      , 0, 'O'},
    { 0, 0, 0, 0}
  };

  setup_logging(argc, argv);

  long requested_seed = 0;
  bool full_obj = false;
  bool face_obj = false;
  unsigned initial_to_beat = 0;
  unsigned lower_bound = 0;
  unsigned improvement_runs = 10;
  unsigned improvement_time = 0;
  unsigned improvement_runs_max = 100000;
  unsigned log_interval = 60;
  int max_time = 0;
  bool obj_in = false;

  while (1) {
    int option_index = 0;
    int r = getopt_long(argc, argv, short_options, long_options, &option_index);

    if (r == -1) break;
    switch (r) {
      case 'h':
        usage(argv[0], 0);
        break;

      case 'S':
        requested_seed = atol(optarg);
        break;

      case 'f':
        full_obj = true;
        break;

      case 'F':
        face_obj = true;
        break;

      case 'b':
        initial_to_beat = atol(optarg);
        break;

      case 'B':
        lower_bound = atol(optarg);
        break;

      case 'I':
        improvement_runs = atol(optarg);
        break;

      case 'i':
        improvement_time = atol(optarg);
        break;

      case 'M':
        improvement_runs_max = atol(optarg);
        break;

      case 'T':
        max_time = atol(optarg);
        break;

      case 'L':
        log_interval = atol(optarg);
        break;

      case 'O':
        obj_in = true;
        break;

      default:
        std::cerr << "Invalid option " << (char)r << std::endl;
        exit(1);
    }
  }

  if (argc - optind > 2) {
    usage(argv[0], 1);
  }

  std::istream *in = &std::cin;
  std::ostream *out = &std::cout;
  std::ifstream filestreamin;
  std::ofstream filestreamout;

  if (argc - optind >= 1) {
    std::string fn(argv[optind]);
    if (fn != "-") {
      filestreamin.open(fn);
      in = &filestreamin;
    }
  }
  if (argc - optind >= 2) {
    std::string fn(argv[optind + 1]);
    if (fn != "-") {
      filestreamout.open(fn);
      out = &filestreamout;
    }
  }

  signal(SIGINT, signalHandler);
  signal(SIGTERM, signalHandler);
  signal(SIGUSR1, signalHandler);
  signal(SIGUSR2, signalHandler);

  auto start_time = std::chrono::system_clock::now();
  auto end_time = start_time + std::chrono::seconds(max_time);
  auto last_info_time = std::chrono::system_clock::now();
  auto info_interval = std::chrono::seconds(log_interval);
  std::chrono::time_point<std::chrono::system_clock> solution_found_at;

  unsigned num_iters = 0;
  unsigned num_iters_since_improved = 0;
  int best_seed = 0;

  bool have_solution = false;
  std::random_device real_rng("/dev/urandom");
  int seed = requested_seed != 0 ? requested_seed : real_rng();
  std::cout << "version: " << GITVERSION << std::endl;
  std::cout << "random_seed: " << seed << std::endl << std::flush;
  random_engine.seed(seed);

  std::unique_ptr<DECL> decl;
  if (obj_in) {
    std::pair<VertexList, InputEdgeSet> p = load_obj(*in);
    decl = std::make_unique<DECL>( std::move(p.first), &p.second );
  } else {
    decl = std::make_unique<DECL>( load_vertices(*in) );
  }
  initial_to_beat = initial_to_beat ? initial_to_beat : decl->get_num_faces();
  unsigned to_beat = initial_to_beat;
  while (1) {
    decl->assert_valid();
    decl->find_convex_decomposition();
    decl->assert_valid();

    ++num_iters;
    ++num_iters_since_improved;

    auto now = std::chrono::system_clock::now();
    unsigned this_num_faces = decl->get_num_faces();

    if (this_num_faces < to_beat) {
      have_solution = true;
      num_iters_since_improved = 0;
      solution_found_at = now;
      to_beat = this_num_faces;
    }

    if (now > last_info_time + info_interval) {
      if (have_solution) {
        LOG(INFO) << "Iter " << num_iters << " overall and " << num_iters_since_improved << "/" << improvement_runs << " since improved; This num faces " << this_num_faces << "; to beat: " << to_beat << "; initial to beat: " << initial_to_beat;
      } else {
        LOG(INFO) << "Iter " << num_iters << " overall; This num faces " << this_num_faces << "; to beat: " << to_beat;
      }
      last_info_time = now;
    }

    if (UNLIKELY(this_num_faces <= lower_bound)) {
      LOG(INFO) << "We hit the lower bound of " << lower_bound;
      std::cout << "exit_reason: lower_bound" << std::endl;
      break;
    } else if (UNLIKELY(have_solution && (num_iters_since_improved >= improvement_runs) && (now > solution_found_at + std::chrono::seconds(improvement_time)))) {
      LOG(INFO) << "We ran for " << num_iters << " overall and " << num_iters_since_improved << "/" << improvement_runs << " since improved.  We did improve on " << initial_to_beat << " faces by " << (initial_to_beat - this_num_faces);
      std::cout << "exit_reason: found-after-improvement_runs" << std::endl;
      std::cout << "improvement_runs: " << improvement_runs << std::endl;
      std::cout << "improvement_time: " << improvement_time << std::endl;
      break;
    } else if (UNLIKELY(!have_solution && num_iters_since_improved >= improvement_runs_max)) {
      LOG(INFO) << "We ran for " << num_iters_since_improved << " without improving on " << initial_to_beat;
      std::cout << "exit_reason: not-found-after-improvement_runs-max" << std::endl;
      break;
    } else if (UNLIKELY(max_time != 0 && now > end_time)) {
      LOG(INFO) << "We ran for max-time of " << max_time << " seconds.  (We did " << num_iters << " iterations total.)";
      std::cout << "exit_reason: timeout" << std::endl;
      break;
    } else if (UNLIKELY(main_loop_interrupted)) {
      std::cout << "exit_reason: interrupt" << std::endl;
      break;
    }
  }

  DBG(DBG_GENERIC) << "Random seed was" << seed;
  std::cout << "num_cvx_areas: " << decl->get_num_faces() << std::endl;
  std::cout << "num_iters: " << num_iters << std::endl;
  auto milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - solution_found_at);
  std::cout << "time_since_found: " << milliseconds.count()/1000. << std::endl;
  milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start_time);
  std::cout << "run_time: " << milliseconds.count()/1000. << std::endl;
  decl->write_obj_segments(full_obj, face_obj, *out);
  return 0;
}
