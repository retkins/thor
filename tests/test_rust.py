import thor 

def test_rust():
    """ Confirm that the Rust backend is loaded 
    """

    y = thor.test(2)
    assert(y == 4.0)