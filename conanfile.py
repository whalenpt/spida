from conan import ConanFile
from conan.tools.cmake import CMakeToolchain, CMakeDeps

class SpidaConan(ConanFile):
    name = "spida"
    version = "0.1"
    settings = "os", "arch", "compiler", "build_type"

    requires = [
        "spdlog/1.17.0",
        "kissfft/131.1.0",
        "boost/1.90.0",
        "nlohmann_json/3.11.3",
    ]

    test_requires = [
        "gtest/1.17.0",
    ]

    default_options = {
        "boost/*:shared": False,
        "gtest/*:shared": False,
        "spdlog/*:shared": False,
        "kissfft/*:shared": False,
        "kissfft/*:datatype": "double",
    }

    def generate(self):
        tc = CMakeToolchain(self)
        tc.cache_variables["CMAKE_EXPORT_COMPILE_COMMANDS"] = True
        tc.generate()

        deps = CMakeDeps(self)
        deps.generate()
