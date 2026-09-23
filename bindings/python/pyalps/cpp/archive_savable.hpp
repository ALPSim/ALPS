// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Declares which bound types can be written with `archive[path] = object`.
//
// pyalps.hdf5's __setitem__ has to decide whether an object knows how to save
// itself or should be picked apart as a number / sequence / mapping. For
// objects defined in Python it duck-types on a bound save() method. That test
// cannot work for extension types: nanobind's bound methods are their own type
// (nb_bound_method), and more importantly having a method *named* save is not
// enough -- alps::alea's MCScalarData / MCVectorData and the pyalea
// observables all have one that takes a file name and an observable name
// instead of an archive, and calling those with an archive is a type error.
//
// So the capability is declared, not guessed, and it is declared at the
// binding site where the signature is known. The archive module then needs no
// knowledge of the types at all -- no includes, no enumeration to keep in sync
// as types are added.
#ifndef PYALPS_ARCHIVE_SAVABLE_HPP
#define PYALPS_ARCHIVE_SAVABLE_HPP
#include <nanobind/nanobind.h>
namespace pyalps {
// Attribute set on the class object. Looked up on instances, so inheritance
// carries it to Python subclasses -- for those, getattr finds the subclass's
// own save() override, which is what should run.
inline constexpr char const * archive_savable_attr = "_alps_archive_savable";
// Mark `cls` as having a save(alps::hdf5::archive &). Call it on every class
// that binds one; the class must also bind "save". Takes a handle so it
// accepts either a nb::class_ or the class looked up on the module, and no
// binding site has to be restructured to name its class object.
inline void mark_archive_savable(nanobind::handle cls) {
    cls.attr(archive_savable_attr) = true;
}
} // namespace pyalps
// Scope: pyalps-internal. This header is not installed, so a downstream
// nanobind module cannot mark its own types -- and does not need to. A class
// exported with ALPS_EXPORT_SIM_TO_PYTHON derives from alps::mcbase and
// inherits the marker from it, and export_sim_to_python binds save() on the
// derived type, so `archive[path] = simulation` reaches the derived save and
// writes the same tree that simulation.save(archive) does. Tutorial 5's
// smoke test asserts exactly that. Ship this header with the SDK only if a
// downstream type ever needs the marker without deriving from mcbase.
#endif // PYALPS_ARCHIVE_SAVABLE_HPP
