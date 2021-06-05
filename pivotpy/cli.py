import pivotpy
import webbrowser
def main():
    print('Pivotpy\n=======')
    print('Version: ',pivotpy.__version__)
    print(pivotpy.__doc__)
    print('Loading Online DOCS...')
    webbrowser.open('https://massgh.github.io/pivotpy/',new=1)