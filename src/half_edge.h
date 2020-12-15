#include "utils.h"
#include <iostream>
#include <math.h>
#include <map>
#include<bits/stdc++.h> 
extern std::map<std::pair<int ,int>,struct halfedge *> dictionary_edges;
extern std::vector<p2t::Point*> points;
typedef struct halfedge{
    struct halfvertex *v;
    struct halfface *f;
    struct halfedge *opposite;
    struct halfedge *next;
    struct halfedge *prev;
} edge;

typedef struct halfface{
    edge *e;
    int triangleType; //2 for terminal, 1 for sleeve and 0 for junction
    int visit;
} face;

typedef struct halfvertex{
    float x;
    float y;
    float z;
    struct halfedge *e;
    int vNum;
} vertex;

static void deleteFaces(std::vector<face*> &faces) {

    int totalFaces = faces.size();

    for(int i = 0; i < totalFaces; ++i)
    {
        delete faces[i]->e->next->next;
        delete faces[i]->e->next;
        delete faces[i]->e;
        delete faces[i];
    }
    faces.clear();
}

static vertex * makeHalfEdgeVertex(float x, float y, float z, int vertexNumber){
    vertex *v = new vertex;
    v->x = x; v->y = y; v->z = z;
    v->e = NULL;  v->vNum = vertexNumber;
    return v;
}

static void makeHalfEdgeFace(int i1,int i2,int i3, const std::vector<vertex *> &vertices, std::vector<face *> &faces){
    // std::cout<<<<std::endl;
    edge *edges[3];
    face *f = new face;
    f->visit = 0;
    f->triangleType = -1;
    edges[0] = new edge;
    edges[1] = new edge;
    edges[2] = new edge;
    f->e = edges[0];

    for (int i =0;i<3;i++){
        edges[i]->f = f;
        edges[i]->opposite = NULL;
        edges[i]->next = edges[(i+1)%3];
        edges[i]->prev = edges[ (i-1)>0 ? i-1 : 2];
    }
    glm::vec3 dir1(vertices[i2]->x-vertices[i1]->x, vertices[i2]->y-vertices[i1]->y, 0.0);
    
    glm::vec3 dir2(vertices[i3]->x-vertices[i1]->x, vertices[i3]->y-vertices[i1]->y, 0.0);

    GLfloat cross_prod_val = glm::cross(dir1, dir2).z;

    if (cross_prod_val<0.0){
        swapInts(i1,i3);
    }
    int indices[] = {i1,i2,i3};
    std::pair<int,int> p;
    for (int i=0;i<3;i++){
        edges[i]->v = vertices.at(indices[i]);
        vertices.at(indices[i])->e = edges[i];

        p.first = std::min(indices[i],indices[(i+1)%3]);
        p.second = std::max(indices[i],indices[(i+1)%3]);
        if (dictionary_edges.find(p) == dictionary_edges.end()){
            dictionary_edges[p] = edges[i];
        }
        else{
            dictionary_edges[p]->opposite = edges[i];
            edges[i]->opposite = dictionary_edges[p];
        }
    }    

    // for (int i =0;i<3)

    faces.push_back(f);
}

static void markTriangles(std::vector<face *>& faces){
    int count;
    edge *e;
    edge *prevE;
    for (auto face:faces){
        count = 0;
        prevE = face->e;
        if (prevE -> opposite == NULL){
            count +=1;
        }
        e = prevE ->next;
        while (e != prevE){
            if (e->opposite == NULL){
                count +=1;
            }
            e = e -> next;
        }
        face->triangleType = count;
    }
}

static void createHalfEdgeBuffers(std::vector<p2t::Point*> &points, std::vector<p2t::Triangle*> &triangles,std::vector<vertex *> &vertices, std::vector<face *> &faces){

    // ;
    std::vector<int> indices;
    float x,y,z;
    float d;
    int indx;
    bool foundIndx;

    for (auto point: points){
        // std::cout<< vertices.size()<<std::endl;
        vertices.push_back(makeHalfEdgeVertex(point->x,point->y,0,vertices.size()));
    }

    for (auto triangle:triangles){
        for (int i =0;i<3;i++){
            x = triangle->GetPoint(i)->x;
            y = triangle->GetPoint(i)->y;
            foundIndx = false;
            indx= 0;
            // std::cout<<"TRIANGGLE"<<std::endl;
            for (auto point: points){
                // std::cout<<x<<", "<<y<<": "<<point->x<<", "<<point->y<<std::endl;
                d = sqrt(pow( x*x - point->x*point->x , 2 )+pow( y*y - point->y*point->y , 2 ));
                // std::cout<< d<<std::endl;
                if (d<=0.00001){
                    // std::cout<<"ok nishant"<<std::endl;
                    foundIndx = true;
                    break;
                }
                indx++;
            }
            indices.push_back(indx);
        }
        makeHalfEdgeFace(indices[0],indices[1],indices[2],vertices,faces);
        indices.clear();
    }
}

static int indexOfVertex(float x, float y, std::vector<vertex *> *vertices){
    int indx=0;
    bool foundIndx=false;
    int centre;
    float d;
    for (auto point: *vertices){
        d = sqrt(pow( x*x  - point->x*point->x , 2 )+pow( y*y  - point->y*point->y , 2 ));
        if (d<=0.00001){
            foundIndx = true;
            break;
        }
        indx++;
    }
    if (!foundIndx){
        centre = vertices->size();
        vertices->push_back(makeHalfEdgeVertex(x,y,0.0,vertices->size()));
    }
    else{
        std::cout<<"--AM I HERE?\n";
        centre = indx;
    }
    std::cout<<"**************************************************"<<vertices->size()<<" "<<centre<<"\n\n";
    return centre;
}

static bool outsideCircle(vertex *v3,float centerX,float centerY, float radius){
    return pow(v3->x-centerX,2)+pow(v3->y-centerY,2)>pow(radius,2);
}

static std::vector<face *> pruneTriangles(std::vector<vertex *> vertices, std::vector<face *> faces){
    edge *prevE;
    edge *e;
    edge *opp;
    int indx;
    float fanPoint[3];
    float d;
    int centre;
    bool outside=false;
    bool foundIndx =false;
    int cnt;
    vertex *v1;
    vertex *v2;
    vertex *v3;
    std::vector<face *> pruned_faces; 
    std::deque<vertex *> uneccVertices;
    for (auto face:faces){
        if (face->triangleType==2){ //terminal
            prevE = face->e; 
            e = face->e;
            uneccVertices.clear();
            outside=false;
            while (!outside){
                std::cout<<"pls"<<std::endl;
                e = e->next;
                cnt =0;
                while (e->opposite==NULL && cnt<=3){
                    e = e->next;
                    std::cout<<"OPPOSITE EDGE: "<<e->opposite<<std::endl;
                    std::cout<<e<<"    "<<prevE<<std::endl;
                    cnt+=1;
                }
                std::cout<<"edge selected:  "<<e<<"  OPPOSITE EDGE: "<<e->opposite<<std::endl;
                face->visit=1;
                e->f->visit=1;
                v1 = e->v;
                v2 = e->next->v;
                v3 = e->next->next->v;
                if (e->opposite->f->triangleType ==0){
                    opp = e ->opposite;
                    fanPoint[0] = (v1->x+v2->x + opp->next->next->v->x)/3;
                    fanPoint[1] = (v1->y+v2->y + opp->next->next->v->y)/3;
                    fanPoint[2] = 0.0;
                    centre = indexOfVertex(fanPoint[0],fanPoint[1],&vertices);
                    if (count(uneccVertices.begin(), uneccVertices.end(), v2) == 0)
                    uneccVertices.push_front(v2);

                    if (count(uneccVertices.begin(), uneccVertices.end(), v1) == 0)
                        uneccVertices.push_back(v1);
                    break;
                }


                if (count(uneccVertices.begin(), uneccVertices.end(), v3) == 0)
                    uneccVertices.push_front(v3);
                
                // vertex * fan_center;
                std::cout<<"here?"<<std::endl;
                float centerX = (v1->x+v2->x)/2.0;
                float centerY = (v1->y+v2->y)/2.0;
                float radius = sqrt(pow(v1->x-v2->x,2) + pow(v1->y-v2->y,2))/2.0;
                for (int i=1;i<=uneccVertices.size()-1;i++){
                    vertex * v = uneccVertices.at(i);
                    if (outsideCircle(v,centerX,centerY,radius)){
                        centre = vertices.size();
                        vertices.push_back(makeHalfEdgeVertex(centerX,centerY,0.0,vertices.size()));
                        outside = true;
                    }
                }
                if (count(uneccVertices.begin(), uneccVertices.end(), v2) == 0)
                    uneccVertices.push_front(v2);

                if (count(uneccVertices.begin(), uneccVertices.end(), v1) == 0)
                    uneccVertices.push_back(v1);
                e = e->opposite;
                std::cout<<"loop agen pls"<<std::endl;
                
            }
            std::cout<<"im out"<<std::endl;
            if (uneccVertices.size() !=0){
                for (int i=0;i<uneccVertices.size()-1;i++){
                    // if (uneccVertices.si)
                    std::cout<<"i value: "<<i<<"   "<<(i+1)<<"  max:  " <<uneccVertices.size()<<std::endl;
                    makeHalfEdgeFace(uneccVertices.at(i)->vNum,uneccVertices.at(i+1)->vNum,centre,vertices, pruned_faces);
                }
            }
        }   
    }

    for (auto face: faces){
        if (face->triangleType==0){
            int centroid;            
            e = face->e;
            for (int i =0;i<3;i++){
                fanPoint[0] = (e->v->x + e->next->v->x)/2.0;
                fanPoint[1] = (e->v->y + e->next->v->y)/2.0;
                fanPoint[2] = 0.0;
                foundIndx = false;
                v1 = e->v;
                v2 = e->next->v;
                v3 = e->next->next->v;
                for (auto point: vertices){
                    d = sqrt(pow( fanPoint[0]*fanPoint[0]  - point->x*point->x , 2 )+pow( fanPoint[1]*fanPoint[1]  - point->y*point->y , 2 ));
                    if (d<=0.00001){
                        foundIndx = true;
                        break;
                    }
                }
                if(e->opposite!=NULL && (e->opposite->f->triangleType==0 || e->opposite->f->visit==0 || foundIndx)){

                    fanPoint[0] = (v1->x + v2->x + v3->x)/3.0;
                    fanPoint[1] = (v1->y + v2->y + v3->y)/3.0;
                    fanPoint[2] = 0.0;
                    centroid = indexOfVertex(fanPoint[0],fanPoint[1],&vertices);
                    centre = indexOfVertex((e->v->x+e->next->v->x)/2.0, (e->v->y+e->next->v->y)/2.0,&vertices);
                    makeHalfEdgeFace(e->v->vNum,centre,centroid,vertices,pruned_faces);
                    makeHalfEdgeFace(centre,e->next->v->vNum,centroid,vertices,pruned_faces);
                }
                e = e->next;
            }
        }

        else if (face->visit==0){
            e = face->e;
            e = e->next;
            cnt =0;
            std::cout<<"pls maut"<<std::endl;
            std::cout<<"Triangle type: "<<face->triangleType<<std::endl;
            std::cout<<"OPPOSITE EDGE: "<<e->opposite<<std::endl;
                    // std::cout<<e<<"    "<<prevE<<std::endl;
            while (e->opposite!=NULL && cnt<=5){
                e = e->next;
                std::cout<<"OPPOSITE EDGE: "<<e->opposite<<std::endl;
                std::cout<<e<<"    "<<prevE<<std::endl;
                cnt +=1;
            }
            face->visit=1;
            e->f->visit=1;
            v1 = e->v;
            v2 = e->next->v;
            v3 = e->next->next->v;
            int center1 = indexOfVertex((v2->x+v3->x)/2.0,(v2->y+v3->y)/2.0,&vertices);
            int center2 = indexOfVertex((v1->x+v3->x)/2.0,(v1->y+v3->y)/2.0,&vertices);
            makeHalfEdgeFace(v1->vNum,v2->vNum,center1,vertices, pruned_faces);
            makeHalfEdgeFace(v1->vNum,center1,center2,vertices, pruned_faces);
            makeHalfEdgeFace(center1,v3->vNum,center2, vertices, pruned_faces);
        }
    }
    std::cout<<"ended"<<std::endl;
    deleteFaces(faces);
    return pruned_faces;
}


