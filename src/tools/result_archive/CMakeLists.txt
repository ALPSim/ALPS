# Historical SQLite-backed XML result indexing/plotting, not HDF5 archive code.
# Register source ownership without restoring an executable or SQLite dependency.
alps_register_module(NAME tools_result_archive
  EXTRA_FILES
    archive.cpp
    archive_index.cpp archive_index.hpp
    archive_node.cpp archive_node.hpp
    archive_plot.cpp archive_plot.hpp
    archive_sqlite.cpp archive_sqlite.hpp
    archive_xml.cpp archive_xml.hpp
  DEPENDS configuration parapack)
