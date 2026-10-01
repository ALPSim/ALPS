# Consume the upstream LP64 DLL package; no ALPS-specific numerical code.
set(VCPKG_BUILD_TYPE release)
# The upstream distribution is one release C-ABI DLL, usable by either MSVC
# configuration. Its imported target intentionally has only a Release location.
set(VCPKG_POLICY_MISMATCHED_NUMBER_OF_BINARIES enabled)
vcpkg_check_linkage(ONLY_DYNAMIC_LIBRARY)
vcpkg_download_distfile(archive
    URLS "https://github.com/OpenMathLib/OpenBLAS/releases/download/v${VERSION}/OpenBLAS-${VERSION}-woa64-dll.zip"
    FILENAME "OpenBLAS-${VERSION}-woa64-dll.zip"
    SHA512 c07ccdb24b7cbc43a18cd105471da2ce98b3b5c7c780d2ce896b450bea0705a91f634855c3c17aa1ae2b6d57319e86831c55a4185b11528326d4165c2dbdb8b2)
vcpkg_extract_source_archive(source ARCHIVE "${archive}")
file(COPY "${source}/bin" "${source}/include" "${source}/lib"
    DESTINATION "${CURRENT_PACKAGES_DIR}")
vcpkg_cmake_config_fixup(PACKAGE_NAME OpenBLAS CONFIG_PATH lib/cmake/OpenBLAS)
vcpkg_fixup_pkgconfig()
vcpkg_download_distfile(license
    URLS "https://raw.githubusercontent.com/OpenMathLib/OpenBLAS/v${VERSION}/LICENSE"
    FILENAME "OpenBLAS-${VERSION}-LICENSE"
    SHA512 e2f8368edd9d35352e7afccc4351aab8774631967c8adc7d8b6f4d1b58a3dddf63492d1bc39e2b2bfe01d6201592428c3f480246cea804bd352ff28f412b1215)
vcpkg_install_copyright(FILE_LIST "${license}")
