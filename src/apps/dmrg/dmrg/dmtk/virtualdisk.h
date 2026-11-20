// ALPS Project: https://alps.comp-phys.org/
// SPDX-License-Identifier: MIT
#ifndef __DMTK_VDISK_H__
#define __DMTK_VDISK_H__

#include <map>
#include <string>

// In-memory replacement for DMRG scratch files: maps a file name to its
// binary contents. Used by System when use_vdisk is enabled.
class VirtualDisk : public std::map<std::string, std::string> {};

#endif // __DMTK_VDISK_H__
