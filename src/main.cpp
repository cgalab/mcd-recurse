#include "io.h"
#include "geom.h"

#include <fstream>
#include <iostream>
#include <getopt.h>

INITIALIZE_EASYLOGGINGPP
unsigned DBG_INDENT_CTR = 0;

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
  FILE *f = err ? stderr : stdout;

  fprintf(f,"Usage: %s [options] <INPUT> <OUTPUT>\n", progname);
  exit(err);
}

int main(int argc, char *argv[]) {
  const char * const short_options = "h";
  const option long_options[] = {
    { "help"        , no_argument      , 0, 'h'},
    { 0, 0, 0, 0}
  };

  setup_logging(argc, argv);

  while (1) {
    int option_index = 0;
    int r = getopt_long(argc, argv, "hc:", long_options, &option_index);

    if (r == -1) break;
    switch (r) {
      case 'h':
        usage(argv[0], 0);
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

  std::unique_ptr<std::vector<Vertex>> vertexlist = load_vertices(*in);
  DECL decl(*vertexlist);
  std::cout << decl << std::endl;

  return 0;
}
