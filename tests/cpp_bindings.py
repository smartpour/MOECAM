from .cffi_interface import ffi, lib

class ToyCppLib:
    def add(self, a, b):
        return lib.add(a, b)

    def print_message(self, message):
        lib.print_message(message.encode('utf-8'))


