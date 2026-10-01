set(VCPKG_POLICY_EMPTY_PACKAGE enabled)
file(INSTALL "${CURRENT_PORT_DIR}/lapack-config.cmake"
    DESTINATION "${CURRENT_PACKAGES_DIR}/share/lapack")
