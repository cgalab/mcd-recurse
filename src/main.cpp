#include "io.h"
#include "geom.h"

#include "gitversion.h"

#include <fstream>
#include <iostream>
#include <chrono>
#include <getopt.h>

INITIALIZE_EASYLOGGINGPP
unsigned DBG_INDENT_CTR = 0;
std::default_random_engine random_engine;

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
    << "    --seed NUM          seed of the RNG" << std::endl
    << "    --full-obj          also print vertex coordinates to .obj file" << std::endl
    << "    --to-beat     NUM   return when a better answer is found" << std::endl
    << "    --lower-bound NUM   return immediately when this is reached" << std::endl
    << "    --min-runs NUM      Do at least NUM runs/attempts at solving this" << std::endl
    << "    --max-time NUM      Do not start a new run after NUM seconds (overrides min-runs)" << std::endl
  ;
  exit(err);
}

int main(int argc, char *argv[]) {
  const char * const short_options = "hS:fb:B:M:T:";
  const option long_options[] = {
    { "help"        , no_argument      , 0, 'h'},
    { "seed"        , required_argument, 0, 'S'},
    { "full-obj"    , no_argument      , 0, 'f'},
    { "to-beat"     , required_argument, 0, 'b'},
    { "lower-bound" , required_argument, 0, 'B'},
    { "min-runs"    , required_argument, 0, 'M'},
    { "max-time"    , required_argument, 0, 'T'},
    { 0, 0, 0, 0}
  };

  setup_logging(argc, argv);

  long requested_seed = 0;
  bool full_obj = false;
  int to_beat = 0;
  int lower_bound = 0;
  int min_runs = 10;
  int max_time = 0;

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

      case 'b':
        to_beat = atol(optarg);
        break;

      case 'B':
        lower_bound = atol(optarg);
        break;

      case 'M':
        min_runs = atol(optarg);
        break;

      case 'T':
        max_time = atol(optarg);
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


  auto start_time = std::chrono::system_clock::now();
  auto end_time = start_time + std::chrono::seconds(max_time);

  int num_iters = 0;
  int best_seed = 0;
  int best_num_faces = to_beat;
  std::stringstream obj_content;

  bool have_solution = false;
  std::random_device real_rng("/dev/urandom");
  int seed = requested_seed != 0 ? requested_seed : real_rng();
  std::cout << "version: " << GITVERSION << std::endl;
  std::cout << "random_seed: " << seed << std::endl << std::flush;
  random_engine.seed(seed);

  DECL decl( load_vertices(*in) );
  while (1) {
    decl.assert_valid();
    decl.find_convex_decomposition();
    decl.assert_valid();

    ++num_iters;

    int this_num_faces = decl.get_num_faces();
    DBG(DBG_GENERIC) << "This num faces " << this_num_faces;

    if (this_num_faces < best_num_faces || best_num_faces == 0) {
      have_solution = true;

      best_num_faces = this_num_faces;
      std::stringstream().swap(obj_content); // clear obj_content
      decl.write_obj_segments(full_obj, obj_content);
    }

    if ( (this_num_faces <= lower_bound)
      || (have_solution && num_iters >= min_runs)
      || (max_time != 0 && std::chrono::system_clock::now() > end_time)
      ) {
      /*
      std::cerr << "a: " << (this_num_faces <= lower_bound) << std::endl;
      std::cerr << "b: " << (have_solution && num_iters >= min_runs) << std::endl;
      std::cerr << "c: " << (max_time != 0 && std::chrono::system_clock::now() > end_time) << std::endl;
      */
      break;
    };
    decl.reset_constraints();
  }

  DBG(DBG_GENERIC) << "Random seed was" << seed;
  if (best_num_faces > 0) {
    std::cout << "num_cvx_areas: " << best_num_faces << std::endl;
    std::cout << "num_iters: " << num_iters << std::endl;
    *out << obj_content.rdbuf();
    return 0;
  } else {
    std::cerr << "No decomposition found." << std::endl;
    return 1;
  }
}
