import os

from conan import ConanFile
from conan.tools.cmake import CMakeToolchain, CMakeDeps, CMake
from conan.tools.files import get


class SpidaConan(ConanFile):
    name = "spida"
    # Kept in sync with CMakeLists.txt's project(SPIDA VERSION ...) and with
    # the vN.N.N git tags release.yml publishes from — source() below relies
    # on this to build the release-tarball download URL.
    version = "0.2.0"
    package_type = "static-library"
    settings = "os", "arch", "compiler", "build_type"

    # Only what's needed to configure+build spida itself (plus its
    # always-vendored NayukiDCT submodule) inside `conan create` — deliberately
    # excludes test/, demos/, worker/, and the other external/ submodules
    # (kissfft, spida's own googletest, the vendored Boost fallback), since
    # those are either unused during packaging (the four Conan-provided deps
    # are always found first) or not part of the installed package. worker/
    # in particular: cmake.install() below only installs the spida library
    # target, never spida-worker, so this recipe has no use for it even
    # though .github/workflows/release.yml's published source tarball
    # includes worker/ (for a different audience — see that file's header
    # comment). That tarball's staged file list intentionally no longer
    # matches this tuple one-for-one; this stays library-scoped on purpose.
    exports_sources = (
        "CMakeLists.txt",
        "cmake/*",
        "include/*",
        "src/*",
        "external/nayukidct/CMakeLists.txt",
        "external/nayukidct/cmake/*",
        "external/nayukidct/include/*",
        "external/nayukidct/src/*",
    )

    test_requires = [
        "gtest/1.17.0",
    ]

    default_options = {
        "boost/*:shared": False,
        "gtest/*:shared": False,
        "spdlog/*:shared": False,
        "kissfft/*:shared": False,
        "kissfft/*:datatype": "double",
        "libiconv/*:shared": False,
    }

    def requirements(self):
        self.requires("spdlog/1.17.0")
        # kissfft's public-header footprint is opaque pointers only
        # (kiss_fft_cfg/kiss_fftr_cfg, forward-declared in transform/fftRVX.h
        # and friends — see those headers) — no transitive_headers needed.
        self.requires("kissfft/131.1.0")
        self.requires("boost/1.90.0")
        # Installed public headers use nlohmann::json by value throughout
        # deeply templated inline code — not an opaque handle, can't be
        # forward-declared away.
        self.requires("nlohmann_json/3.11.3", transitive_headers=True)
        # boost's Boost.Locale iconv codecvt backend pulls libiconv on every
        # non-Linux platform (glibc ships iconv(3) built in, so Linux never
        # needs it) — declared explicitly so it's a visible, direct
        # requirement rather than an invisible transitive one, and so it can
        # be pinned in conan.lock for the macOS CI build.
        if self.settings.os != "Linux":
            self.requires("libiconv/1.17")

    def source(self):
        # exports_sources already populated this when the recipe was
        # exported/created from an actual git checkout of spida (the normal
        # `conan create .` case, verified as part of local dev and CI). Only
        # a bare copy of conanfile.py — no local src/ etc. alongside it, e.g.
        # a consumer without a spida git checkout — falls through here, and
        # fetches the tagged release source tarball instead.
        marker = os.path.join(self.source_folder, "src", "SpidaCVT.cpp")
        if os.path.exists(marker):
            return

        version = str(self.version)
        url = (
            "https://github.com/whalenpt/spida/releases/download/"
            f"v{version}/spida-{version}-src.tar.gz"
        )
        get(self, url, strip_root=True)

    def generate(self):
        tc = CMakeToolchain(self)
        tc.cache_variables["CMAKE_EXPORT_COMPILE_COMMANDS"] = True
        tc.generate()

        deps = CMakeDeps(self)
        deps.generate()

    def build(self):
        # NOTE: deliberately not touching tc.cache_variables in generate()
        # for SPIDA_TEST/SPIDA_DEMOS — that toolchain file is shared with the
        # plain `conan install -of build/<cfg>` + `cmake --preset` flow
        # (CLAUDE.md's mandated build sequence), and a FORCE'd cache_variables
        # entry there would silently clobber the -DSPIDA_TEST=ON/-DSPIDA_DEMOS=ON
        # passed on that later, separate cmake --preset invocation. Both
        # options already default OFF at the top-level CMakeLists.txt, so the
        # package build below needs no overrides at all.
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def package(self):
        cmake = CMake(self)
        cmake.install()

    def package_info(self):
        # spida ships its own hand-authored CMake package config
        # (cmake/SPIDAconfig.cmake.in) that already runs find_dependency()
        # for every one of its own dependencies. Tell Conan not to generate a
        # competing spida-config.cmake via CMakeDeps — just add the installed
        # config's directory to CMAKE_PREFIX_PATH so a consumer's plain
        # find_package(spida) picks up the real one.
        self.cpp_info.set_property("cmake_find_mode", "none")
        self.cpp_info.builddirs.append(os.path.join("lib", "cmake", "spida"))
        # NayukiDCT is a bundled/vendored submodule (external/nayukidct),
        # always built+installed alongside spida rather than pulled in as its
        # own Conan requirement — its config lands in the same package, one
        # directory over, and SPIDAconfig.cmake find_dependency()s it, so it
        # needs to be discoverable via CMAKE_PREFIX_PATH too.
        self.cpp_info.builddirs.append(os.path.join("lib", "cmake", "nayukidct"))
