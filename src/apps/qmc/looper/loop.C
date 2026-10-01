/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1997-2010 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include <looper/version.h>
#ifdef HAVE_PARAPACK_13
# include <alps/parapack/scheduler.h>
#else
# include <alps/parapack/parapack.h>
#endif

int main(int argc, char** argv) { return alps::parapack::start(argc, argv); }

PARAPACK_SET_COPYRIGHT(LOOPER_COPYRIGHT)
PARAPACK_SET_VERSION(LOOPER_VERSION_STRING)
