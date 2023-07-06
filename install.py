import subprocess

def install_libraries():
    libraries = ['pandas', 'numpy', 'matplotlib', 'tkinter']
    
    for library in libraries:
        try:
            __import__(library)
            print(f"{library} is already installed.")
        except ImportError:
            print(f"{library} is not installed. Installing now...")
            subprocess.call(['pip3', 'install', library])
            print(f"{library} has been installed.")
            
install_libraries()
            
            