import os

from conan import ConanFile
from conan.tools.cmake import CMake, cmake_layout
from conan.tools.build import can_run


class SpidaTestPackageConan(ConanFile):
    settings = "os", "arch", "compiler", "build_type"
    generators = "CMakeToolchain", "CMakeDeps"

    def requirements(self):
        # Consumes the spida package `conan create` just built, exactly as an
        # external consumer would — this is what exercises SPIDAconfig.cmake
        # (find_package(spida) + SPIDA::spida) through Conan's dependency
        # graph, not just the plain CMake install tree.
        self.requires(self.tested_reference_str)

    def layout(self):
        cmake_layout(self)

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def test(self):
        if can_run(self):
            cmd = os.path.join(self.cpp.build.bindir, "spida_test_package")
            self.run(cmd, env="conanrun")
