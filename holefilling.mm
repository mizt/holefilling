const bool VERBOSE = false;

#import <Cocoa/Cocoa.h>
#import <simd/simd.h>
#import "holefilling.h"

#import <vcg/space/triangle3.h>
#import <vcg/complex/complex.h>
#import <vcg/complex/algorithms/hole.h>
#import <vcg/complex/algorithms/local_optimization.h>
#import <vcg/complex/algorithms/local_optimization/tri_edge_flip.h>
#import <vcg/complex/algorithms/smooth.h>
#import <vcg/complex/algorithms/refine.h>

bool isSameClassName(id a, NSString *b) { return (a&&[[a className] compare:b]==NSOrderedSame); }
bool isNumber(id a) { return isSameClassName(a,@"__NSCFNumber"); }
bool isBoolean(id a) { return isSameClassName(a,@"__NSCFBoolean"); }

class MyFace;
class MyVertex;
struct MyUsedTypes : public vcg::UsedTypes<vcg::Use<MyVertex>::AsVertexType,vcg::Use<MyFace>::AsFaceType> {};

class MyVertex : public vcg::Vertex<MyUsedTypes, vcg::vertex::Coord3f, vcg::vertex::BitFlags, vcg::vertex::Normal3f, vcg::vertex::Mark, vcg::vertex::Color4b> {};
class MyFace : public vcg::Face<MyUsedTypes, vcg::face::VertexRef, vcg::face::FFAdj, vcg::face::Mark, vcg::face::BitFlags, vcg::face::Normal3f> {};

class MyMesh : public vcg::tri::TriMesh<std::vector<MyVertex>,std::vector<MyFace>> {};

//Delaunay
class MyDelaunayFlip: public vcg::tri::TriEdgeFlip< MyMesh, MyDelaunayFlip > {
public:
    typedef  vcg::tri::TriEdgeFlip< MyMesh,  MyDelaunayFlip > TEF;
    inline MyDelaunayFlip(const TEF::PosType &p, int i, vcg::BaseParameterClass *pp) :TEF(p,i,pp){}
};

bool callback(int percent, const char *str) {
    std::cout << "str: " << str << " " << percent << "%\r";
    return true;
}

template <class MESH>
bool NormalTest(typename vcg::face::Pos<typename MESH::FaceType> pos) {
    //giro intorno al vertice e controllo le normali
    typename MESH::ScalarType thr = 0.0f;
    typename MESH::CoordType NdP = vcg::TriangleNormal<typename MESH::FaceType>(*pos.f);
    typename MESH::CoordType tmp, oop, soglia = typename MESH::CoordType(thr,thr,thr);
    vcg::face::Pos<typename MESH::FaceType> aux=pos;
    do {
        aux.FlipF();
        aux.FlipE();
        oop = Abs(tmp - ::vcg::TriangleNormal<typename MESH::FaceType>(*pos.f));
        if(oop < soglia )return false;
    } while(aux != pos && !aux.IsBorder());
    
    return true;
}

void holefilling(std::vector<simd::float3> *vercites, std::vector<simd::uint3> *faces, NSString *params) {
    
    MyMesh mesh;
    MyMesh::VertexIterator vit = vcg::tri::Allocator<MyMesh>::AddVertices(mesh,vercites->size());
    
    for(int n=0; n<vercites->size(); n++) {
        vit[n].P() = vcg::Point3f(
            (*vercites)[n].x,
            (*vercites)[n].y,
            (*vercites)[n].z
        );
    }
    
    for(int n=0; n<faces->size(); n++) {
        vcg::tri::Allocator<MyMesh>::AddFace(
            mesh,
            &vit[(*faces)[n].x],
            &vit[(*faces)[n].y],
            &vit[(*faces)[n].z]
        );
    }
    
    /*
    " 1) Trivial Ear \n"
    " 2) Minimum weight Ear \n"
    " 3) Selfintersection Ear \n"
    " 4) Minimum weight \n"
    */
    
    int algorithm = 3;
    int holeSize  = 30;
    
    //update the face-face topology
    vcg::tri::UpdateTopology<MyMesh>::FaceFace(mesh);
    vcg::tri::UpdateNormal<MyMesh>::PerVertexPerFace(mesh);
    vcg::tri::UpdateFlags<MyMesh>::FaceBorderFromFF(mesh);
    assert(vcg::tri::Clean<MyMesh>::IsFFAdjacencyConsistent(mesh));

    //compute the average of face area
    float AVG,sumA=0.0f;
    int numA=0,indice;
    indice = mesh.face.size();
    MyMesh::FaceIterator fi;
    for(fi=mesh.face.begin();fi!=mesh.face.end();++fi) {
        sumA+=DoubleArea(*fi)/2;
        numA++;
        for(int ind =0;ind<3;++ind) fi->V(ind)->InitIMark();
    }
    AVG=sumA/numA;
    
    switch(algorithm) {
        case 1: 
            vcg::tri::Hole<MyMesh>::EarCuttingFill<vcg::tri::TrivialEar<MyMesh>>(mesh,holeSize,false); 
            break;
        case 2: 
            vcg::tri::Hole<MyMesh>::EarCuttingFill<vcg::tri::MinimumWeightEar< MyMesh>>(mesh,holeSize,false,callback); 
            break;
        case 3: 
            vcg::tri::Hole<MyMesh>::EarCuttingIntersectionFill<vcg::tri::SelfIntersectionEar< MyMesh>>(mesh,holeSize,false); 
            break;
        case 4:
            vcg::tri::Hole<MyMesh>::MinimumWeightFill(mesh,holeSize,false); 
            vcg::tri::UpdateTopology<MyMesh>::FaceFace(mesh); 
            break;
    }
    
    vcg::tri::UpdateFlags<MyMesh>::FaceBorderFromFF(mesh);
    assert(vcg::tri::Clean<MyMesh>::IsFFAdjacencyConsistent(mesh));
    //printf("\nStart refinig...\n");
    
    /*start refining */
    MyMesh::VertexIterator vi;
    MyMesh::FaceIterator f;
    std::vector<MyMesh::FacePointer> vf;
    f = mesh.face.begin();
    f+=indice;
    for(; f!=mesh.face.end(); ++f) {
        if(!f->IsD()) {
            f->SetS();
        }
    }
    
    std::vector<MyMesh::FacePointer *> FPP;
    std::vector<MyMesh::FacePointer> added;
    std::vector<MyMesh::FacePointer>::iterator vfit;
    int i=1;
    printf("\n");

    for(f=mesh.face.begin(); f!=mesh.face.end(); ++f) {
        if(!(*f).IsD()) {
            if(f->IsS()) {
                f->V(0)->IsW();
                f->V(1)->IsW();
                f->V(2)->IsW();
            }
            else {
                f->V(0)->ClearW();
                f->V(1)->ClearW();
                f->V(2)->ClearW();
            }
        }
    }
    
    vcg::BaseParameterClass pp;
    vcg::LocalOptimization<MyMesh> Fs(mesh,&pp);
    Fs.SetTargetMetric(0.0f);
    Fs.Init<MyDelaunayFlip>();
    Fs.DoOptimization();
    
    do {
        vf.clear();
        f = mesh.face.begin();
        f+=indice;
        for(; f!=mesh.face.end(); ++f) {
            if(f->IsS()) {
                bool test= true;
                for(int ind=0; ind<3; ++ind) f->V(ind)->InitIMark();
                test = (DoubleArea<MyMesh::FaceType>(*f)/2) > AVG;
                if(test) {
                    vf.push_back(&(*f));
                }
            }
        }
        
        //info print
        //printf("\r Refining [%d] - > %d",i,int(vf.size()));
        i++;
        
        FPP.clear();
        added.clear();
        
        for(vfit=vf.begin(); vfit!=vf.end();++vfit) {
            FPP.push_back(&(*vfit));
        }
        int toadd= vf.size();
        MyMesh::FaceIterator f1,f2;
        f2 = vcg::tri::Allocator<MyMesh>::AddFaces(mesh,(toadd*2),FPP);
        MyMesh::VertexIterator vertp = vcg::tri::Allocator<MyMesh>::AddVertices(mesh,toadd);
        std::vector<MyMesh::FacePointer> added;
        added.reserve(toadd);
        vfit=vf.begin();
        
        for(int i=0; i<toadd; ++i,f2++,vertp++) {
            f1=f2;
            f2++;
            vcg::tri::TriSplit<MyMesh,vcg::tri::CenterPointBarycenter<MyMesh> >::Apply(vf[i],&(*f1),&(*f2),&(*vertp),vcg::tri::CenterPointBarycenter<MyMesh>() );
            f1->SetS();
            f2->SetS();
            for(int itr=0; itr<3; itr++) {
                f1->V(itr)->SetW();
                f2->V(itr)->SetW();
            }
            added.push_back(&(*f1));
            added.push_back(&(*f2));
        }
        
        vcg::BaseParameterClass pp;
        vcg::LocalOptimization<MyMesh> FlippingSession(mesh,&pp);
        FlippingSession.SetTargetMetric(0.0f);
        FlippingSession.Init<MyDelaunayFlip >();
        FlippingSession.DoOptimization();
            
    } while(!vf.empty());
    
    vcg::LocalOptimization<MyMesh> Fiss(mesh,&pp);
    Fiss.SetTargetMetric(0.0f);
    Fiss.Init<MyDelaunayFlip>();
    Fiss.DoOptimization();
    
    int UBIT = MyMesh::VertexType::NewBitFlag();
    f=mesh.face.begin();
    f+=indice;
    for(; f!=mesh.face.end(); ++f) {
        if(f->IsS()) {
            for(int ind=0; ind<3; ++ind) {
                if(NormalTest<MyMesh>(vcg::face::Pos<MyMesh::FaceType>(&(*f),ind))) {
                    f->V(ind)->SetUserBit(UBIT);
                }
            }
            f->ClearS();
        }
    }
    
    for(vi=mesh.vert.begin(); vi!=mesh.vert.end(); ++vi) {
        if(!(*vi).IsD()) {
            if(vi->IsUserBit(UBIT)) {
                (*vi).SetS();
                vi->ClearUserBit(UBIT);
            }
        }
    }
    
    vcg::tri::Smooth<MyMesh>::VertexCoordLaplacian(mesh,1,true);
    
    vercites->clear();
    faces->clear();
    
    unsigned int num = 0;
    std::vector<int> indices(mesh.vert.size());
    for(unsigned int n=0; n<mesh.vert.size(); n++) {
        if(!mesh.vert[n].IsD()) {
            indices[n]=num++;
            vercites->push_back(simd::float3{
                mesh.vert[n].P()[0],
                mesh.vert[n].P()[1],
                mesh.vert[n].P()[2]
            });
        }
    }
    
    for(unsigned int n=0; n<mesh.face.size(); n++) {
        if(!mesh.face[n].IsD()) {
            if(mesh.face[n].VN()==3) {
                faces->push_back(simd::uint3{
                    (unsigned int)indices[vcg::tri::Index(mesh,mesh.face[n].V(0))],
                    (unsigned int)indices[vcg::tri::Index(mesh,mesh.face[n].V(1))],
                    (unsigned int)indices[vcg::tri::Index(mesh,mesh.face[n].V(2))]
                });
            }
        }
    }
}