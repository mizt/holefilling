#import <Cocoa/Cocoa.h>
#import <simd/simd.h>
#import <vector>

namespace FileManager {
  bool exists(NSString *path) {
    static const NSFileManager *fileManager = [NSFileManager defaultManager];
    NSError *err = nil;
    [fileManager attributesOfItemAtPath:path error:&err];
    return (err)?false:true;
  }
}

typedef void (*FILTER)(std::vector<simd::float3> *, std::vector<simd::uint3> *, NSString *);

int main(int argc, char *argv[]) {
  
  NSString *obj = @"./mesh.obj";
  
  NSString *pluginPath = @"./holefilling.plugin";
  CFStringRef pluginFunctionName = (__bridge CFStringRef)(@"holefilling");
  
  NSString *pluginParams = [NSString stringWithFormat:@"%s",R"({})"];
  
  if(FileManager::exists(obj)&&FileManager::exists(pluginPath)) {
    
    NSString *src = [NSString stringWithContentsOfFile:obj encoding:NSUTF8StringEncoding error:nil];
    NSArray *lines = [src componentsSeparatedByCharactersInSet:[NSCharacterSet newlineCharacterSet]];
    NSCharacterSet *whitespaces = [NSCharacterSet whitespaceCharacterSet];
    
    std::vector<simd::float3> v;
    std::vector<simd::uint3> f;
    
    for(int k=0; k<lines.count; k++) {
      NSArray *arr = [lines[k] componentsSeparatedByCharactersInSet:whitespaces];
      if([arr count]>0) {
        if([arr[0] isEqualToString:@"v"]&&[arr count]>=4) {
          v.push_back(simd::float3{
            [arr[1] floatValue],
            [arr[2] floatValue],
            [arr[3] floatValue],
          });
        }
        else if([arr[0] isEqualToString:@"f"]&&[arr count]==4) {
          f.push_back(simd::uint3{
            (unsigned int)([arr[1] intValue]-1),
            (unsigned int)([arr[2] intValue]-1),
            (unsigned int)([arr[3] intValue]-1)
          });
        }
      }
    }
    
    if(FileManager::exists(pluginPath)) {
      //NSLog(@"%d",pluginPass);
      CFBundleRef bundle = CFBundleCreate(kCFAllocatorDefault,(CFURLRef)[NSURL fileURLWithPath:pluginPath]);
      if(bundle) {
        FILTER filter = (FILTER)CFBundleGetFunctionPointerForName(bundle,pluginFunctionName);
        if(filter) {
          filter(&v,&f,pluginParams);
          if(bundle) {
            filter = nullptr;
            CFBundleUnloadExecutable(bundle);
            CFRelease(bundle);
          }
        }
      }
    }
    
    NSMutableString *obj = [NSMutableString stringWithString:@""];
    for(unsigned int n=0; n<v.size(); n++) {
      [obj appendString:[NSString stringWithFormat:@"v %0.7f %0.7f %0.7f\n",v[n].x,v[n].y,v[n].z]];
    }
    for(unsigned int n=0; n<f.size(); n++) {
      [obj appendString:[NSString stringWithFormat:@"f %d %d %d\n",f[n].x+1,f[n].y+1,f[n].z+1]];
    }
    [obj writeToFile:@"./dst.obj" atomically:YES encoding:NSUTF8StringEncoding error:nil];
    
  }
  
  return 0;
}