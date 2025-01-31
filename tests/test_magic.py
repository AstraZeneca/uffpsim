
import uffpsim
import os

def test_magic_number():
    if os.path.isfile("pytest-magic.h5"):
        os.remove("pytest-magic.h5")
    fpstore = uffpsim.uffpsimLib.FingerprintStore('pytest-magic.h5', mode="w")
    assert fpstore.set_magic_number()
    assert fpstore.magic_number_exists()
    if os.path.isfile("pytest-magic.h5"):
        os.remove("pytest-magic.h5")
