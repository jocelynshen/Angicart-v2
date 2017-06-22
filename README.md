# Angicart++

### Installation
*Download XCode*
>xcode-select --install
*Install wxmac*
>brew install wxmac
*Install homebrew*
> /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
*Check for libraries*
> /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

### Using Angicart++
You should now be able to look at the code in the source folder, and import example image files and output data.

1. Make sure the path is right for the folders you want to analyze. Path line is specified near the end of *minSurfTests.xcodeproj* and reads:
  <pre><code>string pathToImages("/Users/vsavage/Desktop/minimalSurfaces_share");</pre></code>
1. Look for the line that reads:
  string imageDir(pathToImages + "/FITC-MCA0_N12_NIH001_s1_4"); int imStart(0), imEnd(66); double voxdims[] = {1.648, 1.648, 0.492};  thresh = 0.22; string lengthUnit("um");
1. Where it says, “/FITC-MCA0_N12_NIH001_s1_4” is where you can replace the folder name for whatever images you want to analyze (ie. “small_easy_test”)
1. Voxel dimensions must be where it reads:
  voxdims[] = {1.648, 1.648, 0.492}
Note: voxdims[] = {x, y, z} where x, y, and z are the 3D dimensions for an image (in units of microns)
1. You can now run *minSurfTests.xcodeproj* by clicking the play button at the top left corner of the XCode window, or going to the "Product" menu and selecting "Run", or simultaneously pressing "Command" + "R"
The movie will run, and the output file for vessel measurements

_Note:_
There is currently no gui for choosing file names or even example files such as small_easy_test of FITC…. folders.
