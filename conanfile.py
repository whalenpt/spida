from conan import ConanFile
from conan.tools.cmake import CMakeToolchain, CMakeDeps, cmake_layout


class SpidaConan(ConanFile):
    name = "spida"
    version = "0.1"
    settings = "os", "arch", "compiler", "build_type"

    # -----------------------------
    # Dependencies (prefer pinned or ~ ranges for reproducibility)
    # -----------------------------
    requires = [
        "spdlog/1.17.0",
        "kissfft/131.1.0",
        "boost/1.90.0",
    ]
    test_requires = [
        "gtest/1.17.0",
    ]

    # -----------------------------
    # Build options (important for consistency)
    # -----------------------------
    default_options = {
        "boost/*:shared": False,
        "gtest/*:shared": False,
        "spdlog/*:shared": False,
        "kissfft/*:shared": False,
    }

    def layout(self):
        cmake_layout(self)

    def generate(self):
        tc = CMakeToolchain(self)
        tc.cache_variables["CMAKE_EXPORT_COMPILE_COMMANDS"] = True
        tc.generate()
        deps = CMakeDeps(self)
        deps.generate()
