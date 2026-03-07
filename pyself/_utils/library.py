import ctypes.util
import os


def find_library_full_path(library_name):
    """
    Finds the full path of a library using ctypes.util.find_library.

    Args:
        library_name: The name of the library to find.

    Returns:
        The full path of the library, or None if not found.
    """
    library_path = ctypes.util.find_library(library_name)
    print(f"Library path for {library_name}: {library_path}")
    if library_path:
        if os.path.isabs(library_path):
            return library_path
        else:
            # On Linux, find_library often returns just the filename, so we search in common library paths
            for path in ["/lib", "/usr/lib", "/usr/local/lib"]:
                full_path = os.path.join(path, library_path)
                if os.path.exists(full_path):
                    return full_path
            # If not found in standard paths, try searching in the directories in LD_LIBRARY_PATH
            ld_library_path = os.environ.get("LD_LIBRARY_PATH")
            print(f"LD_LIBRARY_PATH: {ld_library_path}")
            if ld_library_path:
                for path in ld_library_path.split(":"):
                    full_path = os.path.join(path, library_path)
                    print(
                        f"Checking {full_path} for library {library_name} : {os.path.exists(full_path)}"
                    )
                    if os.path.exists(full_path):
                        return full_path
            # If still not found, return None
            return None
    return None
