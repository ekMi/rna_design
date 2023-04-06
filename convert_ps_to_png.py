from PIL import Image


def convert_to_png(path):
    img = Image.open(path)
    img.save("img.png")

if __name__ == "__main__":
    path = 'test/test_files/test_ss.ps'
    convert_to_png(path)
