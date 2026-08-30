from conan import ConanFile
from conan.tools.cmake import CMakeToolchain, CMakeDeps

class SpidaConan(ConanFile):
    name = "spida"
    version = "0.1"
    settings = "os", "arch", "compiler", "build_type"

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
        self.requires("kissfft/131.1.0")
        self.requires("boost/1.90.0")
        self.requires("nlohmann_json/3.11.3")
        # boost's Boost.Locale iconv codecvt backend pulls libiconv on every
        # non-Linux platform (glibc ships iconv(3) built in, so Linux never
        # needs it) — declared explicitly so it's a visible, direct
        # requirement rather than an invisible transitive one, and so it can
        # be pinned in conan.lock for the macOS CI build.
        if self.settings.os != "Linux":
            self.requires("libiconv/1.17")

    def generate(self):
        tc = CMakeToolchain(self)
        tc.cache_variables["CMAKE_EXPORT_COMPILE_COMMANDS"] = True
        tc.generate()

        deps = CMakeDeps(self)
        deps.generate()
