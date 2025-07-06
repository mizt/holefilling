cd "$(dirname "$0")"
cd ./

set -eu

PLUGIN_NAME="holefilling"
PLUGIN_PATH="./${PLUGIN_NAME}.plugin"

echo ${PLUGIN_PATH}

mkdir -p ${PLUGIN_PATH}/Contents/MacOS
clang++ -std=c++20 -Wc++20-extensions -bundle -fobjc-arc -O3 -I ./ -I./eigen/3.4.0_1/include/eigen3 -framework Cocoa ./${PLUGIN_NAME}.mm \
-o ${PLUGIN_PATH}/Contents/MacOS/${PLUGIN_NAME}
cp ./Info.plist ${PLUGIN_PATH}/Contents/

# codesign --force --options runtime --deep --entitlements "./entitlements.plist" --sign "Developer ID Application" --timestamp --verbose ${PLUGIN}.plugin

echo "** BUILD SUCCEEDED **"
