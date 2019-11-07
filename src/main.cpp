#include "io.h"
#include "geom.h"

#include <fstream>
#include <iostream>
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
    << "    --seed SEED   seed of the RNG" << std::endl
    << "    --full-obj    also print vertex coordinates to .obj file" << std::endl
  ;
  exit(err);
}

int main(int argc, char *argv[]) {
  const char * const short_options = "hS:";
  const option long_options[] = {
    { "help"        , no_argument      , 0, 'h'},
    { "seed"        , required_argument, 0, 'S'},
    { "full-obj"    , no_argument      , 0, 'f'},
    { 0, 0, 0, 0}
  };

  setup_logging(argc, argv);

  long seed = std::random_device("/dev/urandom")();
  bool full_obj = false;

  while (1) {
    int option_index = 0;
    int r = getopt_long(argc, argv, "hc:f", long_options, &option_index);

    if (r == -1) break;
    switch (r) {
      case 'h':
        usage(argv[0], 0);
        break;

      case 'S':
        seed = atol(optarg);
        break;

      case 'f':
        full_obj = true;
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

  random_engine.seed(seed);
  std::cout << "random_seed: " << seed << std::endl << std::flush;

  std::unique_ptr<std::vector<Vertex>> vertexlist = load_vertices(*in);
  DECL decl(*vertexlist);
  decl.assert_valid();
  decl.unconstrain_all();
  decl.assert_valid();
  std::cout << "num_cvx_areas: " << decl.get_num_faces() << std::endl;
  decl.write_obj_segments(full_obj ? &*vertexlist : NULL, *out);

  return 0;
}
