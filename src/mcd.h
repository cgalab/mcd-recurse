/* CG:SHOP 2020: Minimum Convex Decomposition -- Recursor Tool
*
*  Copyright 2019, 2020 Peter Palfraader
*
*  This program is free software: you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation, either version 3 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#pragma once

#include "config.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "easyloggingpp/src/easylogging++.h"
#pragma GCC diagnostic pop

#include <assert.h>
#include <random>

extern std::default_random_engine random_engine;
extern bool main_loop_interrupted;

/** Expands to a syntactically valid empty statement.  */
#define STMT_NIL (void)0
#define STMT_BEGIN do {
#define STMT_END } while (0)


/* count how many steps we needed to find an edge incident to a higher-degree vertex */
/* #define COUNT_SEARCHING_FOR_THE_RIGHT_EDGE */

#define DBG_GENERIC                       ( 1u << 0 )
#define DBG_ITERATORS                     ( 1u << 1 )
#define DBG_SHOOTHOLE                     ( 1u << 2 )
#define DBG_SHOOTHOLE2                    ( 1u << 3 )
#define DBG_DECOMPOSITION_LOOP            ( 1u << 4 )
#define DBG_UNCONSTRAIN                   ( 1u << 5 )
#define DBG_FLIP                          ( 1u << 6 )
// #define DBG_STATESAVE                     ( 1u << x )
// #define DBG_STATERESTORE                  ( 1u << x )

#define DEBUG_MASK (                  \
                DBG_GENERIC                       | \
                DBG_SHOOTHOLE                     | \
                0 )
/*
                DBG_FLIP                          | \
                DBG_UNCONSTRAIN                   | \
                DBG_SHOOTHOLE2                    | \
                DBG_DECOMPOSITION_LOOP            | \

                DBG_ITERATORS                     | \

                DBG_STATESAVE                     | \
                DBG_STATERESTORE                  | \
*/

#ifndef NDEBUG
  #define DEBUG_STMT(x) x
#else
  #define DEBUG_STMT(x) STMT_NIL
#endif


#ifdef DEBUG_OUTPUT
  extern unsigned DBG_INDENT_CTR;
  inline void DBG_INDENT_INC() { ++DBG_INDENT_CTR; }
  inline void DBG_INDENT_DEC() { assert(DBG_INDENT_CTR > 0); --DBG_INDENT_CTR; }
  inline std::string DBG_INDENT() {
    std::ostringstream oss;
    for (unsigned i=0; i<DBG_INDENT_CTR; ++i) {
      oss << (((i+1)%2 == 0)  ? u8"· " : "  ");
    }
    return oss.str();
  }

  #ifdef DEBUG_OUTPUT_WITH_FILES
    #define __FILENAME__ (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__)
    #define LOG_EXTRA_INFO << __FILENAME__ << ":" << __LINE__ << " "
  #else
    #define LOG_EXTRA_INFO
  #endif

  #define DBG_INDENT_LEVEL_STORE unsigned __dbg_indent_ctr = DBG_INDENT_CTR
  #define DBG_INDENT_LEVEL_CHECK assert(__dbg_indent_ctr == DBG_INDENT_CTR)

  #define DBG(FACILITY) if (DEBUG_MASK & (FACILITY)) LOG(DEBUG) << DBG_INDENT()  LOG_EXTRA_INFO << __FUNCTION__ << "(): "
#else
  inline void DBG_INDENT_INC() {}
  inline void DBG_INDENT_DEC() {}

  #define DBG_INDENT_LEVEL_STORE STMT_NIL
  #define DBG_INDENT_LEVEL_CHECK STMT_NIL

  #define DBG(x) do { (void)(x); } while (0); if (0) LOG(DEBUG)
#endif

#define DBG_FUNC_BEGIN(FACILITY) do { DBG_INDENT_INC(); DBG(FACILITY) << "{"; } while (0)
#define DBG_FUNC_END(FACILITY)   do { DBG(FACILITY) << "}"  ; DBG_INDENT_DEC(); } while (0)
