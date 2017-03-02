#include "glCanvas.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <glm/gtx/vector_angle.hpp>
#include <iterator>
#include <set>
#include "mesh.h"
#include "edge.h"
#include "vertex.h"
#include "triangle.h"

//A custom function to calculate average normal for gouraud shading
glm::vec3 Mesh::average_normal(Edge *e){

    glm::vec3 result;

    //Collect triangles connected to the edge's end point
    //except 2 triangles adjacent to the edge
    std::vector<Triangle*> triangles = extract_triangles(e);

    //collect normal associate to edge
    glm::vec3 main_normal = ComputeNormal(e->getTriangle()->operator[](0)->getPos(),
                                          e->getTriangle()->operator[](1)->getPos(),
                                          e->getTriangle()->operator[](2)->getPos());

    //Collect triangle adjacent to the edge
    if (e->getOpposite()!=NULL){
        triangles.push_back(e->getOpposite()->getTriangle());
    }

    //assign first normal

    result = main_normal;

    //Computer the rest normals
    for (int i = 0; i < triangles.size(); ++i) {
        glm::vec3 normal;
        normal = ComputeNormal(triangles[i]->operator[](0)->getPos(),
                                triangles[i]->operator[](1)->getPos(),
                                triangles[i]->operator[](2)->getPos());

        result += normal;
        /*
        if(glm::length(main_normal-normal) < 1.117){

        }
         */

    }

    //return average normal
    return glm::normalize(result);

}

// =================================================================
// SUBDIVISION
// =================================================================


void Mesh::LoopSubdivision() {
    printf("Subdivide the mesh!\n");

    //variable structure for a single triangle
    //

                         Vertex *top;
    //                     [0]
    //                    /   \
    //                   /     \

            Vertex *med_1;     Vertex *med_2;

    //                /           \
    //               /             \

    Vertex *bottom_1;Vertex *bottom_2;Vertex *bottom_3;
    //     [1]                                 [2]


    //Temporary store all old triangle pointers
    std::vector<Triangle *> triangle_temp;
    for (triangleshashtype::iterator iter = triangles.begin(); iter != triangles.end(); iter++) {
        triangle_temp.push_back(iter->second);
    }

    //Make a copy of old vertices
    std::vector<Vertex*> vertices_back_up;


    //Store connection between new and old vertices
    std::unordered_map<Vertex*, Vertex*> vertex_map;
    for (int k = 0; k < vertices.size(); ++k) {
        Vertex *temp = new Vertex(vertices[k]->getIndex(),vertices[k]->getPos());
        vertex_map[ vertices[k] ]=temp;
        vertices_back_up.push_back(temp);
    }

    //add new triangles
    for (int i = 0; i < triangle_temp.size(); ++i) {
        //Collect odd vertices
        top = triangle_temp[i]->operator[](0);
        bottom_1 = triangle_temp[i]->operator[](1);
        bottom_3 = triangle_temp[i]->operator[](2);

        //store 3 creases between odd vertices
        float top_to_bottom1 = 0;
        float bottom1_to_bottom3 = 0;
        float bottom3_to_top = 0;

        if (edges[std::make_pair(top,bottom_1)]->getCrease() > 0){
            top_to_bottom1 = edges[std::make_pair(top,bottom_1)]->getCrease();
        }
        if (edges[std::make_pair(bottom_1,bottom_3)]->getCrease() > 0){
            bottom1_to_bottom3 = edges[std::make_pair(bottom_1,bottom_3)]->getCrease();
        }
        if (edges[std::make_pair(bottom_3,top)]->getCrease() > 0){
            bottom3_to_top = edges[std::make_pair(bottom_3,top)]->getCrease();
        }

        //adjust med_1 accoding to Hoppe's mask
        std::tuple<double, double, double> med_1_t = adjust_odd_vertex(edges[std::make_pair(top, bottom_1)], vertex_map);
        glm::vec3 med1_new = glm::vec3(std::get<0>(med_1_t), std::get<1>(med_1_t), std::get<2>(med_1_t));

        //adjust med_2 accoding to Hoppe's mask
        std::tuple<double, double, double> med_2_t = adjust_odd_vertex(edges[std::make_pair(bottom_3, top)], vertex_map);
        glm::vec3 med2_new = glm::vec3(std::get<0>(med_2_t), std::get<1>(med_2_t), std::get<2>(med_2_t));

        //adjust bottom_2 accoding to Hoppe's mask
        std::tuple<double, double, double> bottom_2_t = adjust_odd_vertex(edges[std::make_pair(bottom_1, bottom_3)], vertex_map);
        glm::vec3 bottom2_new = glm::vec3(std::get<0>(bottom_2_t), std::get<1>(bottom_2_t), std::get<2>(bottom_2_t));


        //adjust top accoding to Hoppe's mask
        std::tuple<double, double, double> top_t = adjust_even_vertex(edges[std::make_pair(bottom_3, top)] , vertex_map);
        glm::vec3 top_new = glm::vec3(std::get<0>(top_t), std::get<1>(top_t), std::get<2>(top_t));
        top->setPos(top_new);

        //adjust bottom_1 accoding to Hoppe's mask
        std::tuple<double, double, double> bottom_1_t = adjust_even_vertex(edges[std::make_pair(top, bottom_1)] , vertex_map);
        glm::vec3 bottom_1_new = glm::vec3(std::get<0>(bottom_1_t), std::get<1>(bottom_1_t), std::get<2>(bottom_1_t));
        bottom_1->setPos(bottom_1_new);

        //adjust bottom_2 accoding to Hoppe's mask
        std::tuple<double, double, double> bottom_3_t = adjust_even_vertex(edges[std::make_pair(bottom_1, bottom_3)] , vertex_map);
        glm::vec3 bottom_3_new = glm::vec3(std::get<0>(bottom_3_t), std::get<1>(bottom_3_t), std::get<2>(bottom_3_t));
        bottom_3->setPos(bottom_3_new);

        //check if any of med_1 med_2 or bottom_2 already exits,
        //if so redirect the pointer to the exiting object
        //else add the new vertex
        bool med_1_exist = false;
        bool med_2_exist = false;
        bool bottom_2_exist = false;

        //check if they already exits
        for (int j = 0; j < vertices.size(); ++j) {
            if (vertices[j]->getPos() == med1_new) {
                med_1_exist = true;
                med_1 = vertices[j];
            }
            if (vertices[j]->getPos() == med2_new) {
                med_2_exist = true;
                med_2 = vertices[j];
            }
            if (vertices[j]->getPos() == bottom2_new) {
                bottom_2_exist = true;
                bottom_2 = vertices[j];
            }
        }
        if (med_1_exist == false) {
            med_1 = addVertex(med1_new);
        }
        if (med_2_exist == false) {
            med_2 = addVertex(med2_new);
        }
        if (bottom_2_exist == false) {
            bottom_2 = addVertex(bottom2_new);
        }


        //add 4 new small triangles
        addTriangle(top, med_1, med_2);
        addTriangle(med_1, bottom_1, bottom_2);
        addTriangle(med_2, bottom_2, bottom_3);
        addTriangle(med_2, med_1, bottom_2);

        //Since we want the crease property exit in next generation
        //we pass them to the next generation to revelant edges
        if (top_to_bottom1!=0){
            edges[std::make_pair(top,med_1)]->setCrease(top_to_bottom1);
            edges[std::make_pair(med_1,bottom_1)]->setCrease(top_to_bottom1);
        }
        if (bottom1_to_bottom3!=0){
            edges[std::make_pair(bottom_1,bottom_2)]->setCrease(bottom1_to_bottom3);
            edges[std::make_pair(bottom_2,bottom_3)]->setCrease(bottom1_to_bottom3);
        }
        if (bottom3_to_top!=0){
            edges[std::make_pair(bottom_3,med_2)]->setCrease(bottom3_to_top);
            edges[std::make_pair(med_2,top)]->setCrease(bottom3_to_top);
        }

    }



    //remove all triangles
    for (int j = 0; j < triangle_temp.size(); ++j) {
        removeTriangle(triangle_temp[j]);

    }

    //remove vertecies backup
    for (int l = 0; l < vertices_back_up.size(); ++l) {
        delete vertices_back_up[l];
    }


}

std::tuple<double ,double ,double> Mesh::adjust_odd_vertex(Edge *e , std::unordered_map<Vertex*,Vertex*> convertion){
    //This function is using Hoppe's mask to adjust odd vertices
    double x = 0;
    double y = 0;
    double z = 0;
    double tbe = 3.0f/8.0f;
    double obe = 1.0f/8.0f;
    double obt = 0.5f;

    if (e->getOpposite()!=NULL && e->getCrease()==0){
        //if it is a regular interior vertex
        x += tbe * convertion[e->getStartVertex()]->x();
        y += tbe * convertion[e->getStartVertex()]->y();
        z += tbe * convertion[e->getStartVertex()]->z();

        x += tbe * convertion[e->getEndVertex()]->x();
        y += tbe * convertion[e->getEndVertex()]->y();
        z += tbe * convertion[e->getEndVertex()]->z();

        x += obe * convertion[e->getNext()->getEndVertex()]->x();
        y += obe * convertion[e->getNext()->getEndVertex()]->y();
        z += obe * convertion[e->getNext()->getEndVertex()]->z();

        x += obe * convertion[e->getOpposite()->getNext()->getEndVertex()]->x();
        y += obe * convertion[e->getOpposite()->getNext()->getEndVertex()]->y();
        z += obe * convertion[e->getOpposite()->getNext()->getEndVertex()]->z();
    }
    else{

        //if it is a bpundary or crease
        x += obt * convertion[e->getStartVertex()]->x();
        y += obt * convertion[e->getStartVertex()]->y();
        z += obt * convertion[e->getStartVertex()]->z();

        x += obt * convertion[e->getEndVertex()]->x();
        y += obt * convertion[e->getEndVertex()]->y();
        z += obt * convertion[e->getEndVertex()]->z();
    }

    return std::make_tuple(x,y,z);
}

std::tuple<double ,double ,double> Mesh::adjust_even_vertex(Edge *e, std::unordered_map<Vertex*,Vertex*> convertion){
    //This function is using Hoppe's mask to adjust odd vertices
    double x = 0;
    double y = 0;
    double z = 0;
    double tbf = 3.0f/4.0f;
    double obe = 1.0f/8.0f;
    double tbe = 3.0f/8.0f;
    double tbs = 3.0f/16.0f;

    //adjust the vertex at the end of the edge
    double valence=1;
    int type = 0; //0 for regular, 1 for boundary
    Edge *boundary_edge1; //if exists denotes from other to edge's end
    Edge *boundary_edge2; //if exists denotes from edge's end to other vertex
    Edge *edge_working_on = e;
    std::vector<Edge*> crease_edges; //a vector of all crease edges (do not include boundry)
    std::vector<Vertex*> around_edges; //a vector of all around edges (do not include crease)

    //identity current edge first
    if(e->getCrease()>0){
        crease_edges.push_back(e);
    }
    around_edges.push_back(e->getStartVertex());


    //trying to identfy if all edges around the vertex are in a circle
    while(edge_working_on->getOpposite()!=NULL){
        edge_working_on = edge_working_on->getOpposite();
        edge_working_on = edge_working_on->getNext()->getNext();
        valence++;
        around_edges.push_back(edge_working_on->getStartVertex());
        if (edge_working_on == e){
            type=0;
            break;
        }
        else{
            if ((edge_working_on->getCrease()>0) && (edge_working_on->getOpposite()!=NULL)){
                crease_edges.push_back(edge_working_on);
            }
        }
    }

    //if not set type to boundary
    if(edge_working_on->getOpposite()==NULL){
        type=1;
        boundary_edge1 = edge_working_on;
    }

    edge_working_on =e;

    if(type==1){
        //if boundary, try to find another edge
        while(edge_working_on->getNext()->getOpposite()!=NULL){
            edge_working_on = edge_working_on->getNext()->getOpposite();
            around_edges.push_back(edge_working_on->getStartVertex());
            valence++;
        }

        edge_working_on = edge_working_on->getNext();
        around_edges.push_back(edge_working_on->getEndVertex());
        boundary_edge2 = edge_working_on;
        valence++;
    }


    //start applying Hoppe's mask
    if (type==0 && crease_edges.size() < 2){
        //if not boundary and no more than two crease edge connected
        double beta;
        if (valence>3.0){
            beta = tbe/valence;
        }
        else{
            beta = tbs;
        }
        for (int i = 0; i < around_edges.size(); ++i) {
            x += convertion[around_edges[i]]->x() * beta;
            y += convertion[around_edges[i]]->y() * beta;
            z += convertion[around_edges[i]]->z() * beta;
        }
        x += (1.0 - valence * beta) * convertion[e->getEndVertex()]->x();
        y += (1.0 - valence * beta) * convertion[e->getEndVertex()]->y();
        z += (1.0 - valence * beta) * convertion[e->getEndVertex()]->z();
    }
    else if (type==0 && crease_edges.size()==2){
        //if not boundary and exactly two crease edges connected
        x += convertion[crease_edges[0]->getStartVertex()]->x() * obe;
        y += convertion[crease_edges[0]->getStartVertex()]->y() * obe;
        z += convertion[crease_edges[0]->getStartVertex()]->z() * obe;

        x += convertion[crease_edges[1]->getStartVertex()]->x() * obe;
        y += convertion[crease_edges[1]->getStartVertex()]->y() * obe;
        z += convertion[crease_edges[1]->getStartVertex()]->z() * obe;

        x += convertion[e->getEndVertex()]->x() * tbf;
        y += convertion[e->getEndVertex()]->y() * tbf;
        z += convertion[e->getEndVertex()]->z() * tbf;

    }
    else if (type==1 && crease_edges.size()==0){
        //if boundary and no crease edges connected
        x += convertion[boundary_edge1->getStartVertex()]->x() * obe;
        y += convertion[boundary_edge1->getStartVertex()]->y() * obe;
        z += convertion[boundary_edge1->getStartVertex()]->z() * obe;

        x += convertion[boundary_edge2->getEndVertex()]->x() * obe;
        y += convertion[boundary_edge2->getEndVertex()]->y() * obe;
        z += convertion[boundary_edge2->getEndVertex()]->z() * obe;

        x += convertion[e->getEndVertex()]->x() * tbf;
        y += convertion[e->getEndVertex()]->y() * tbf;
        z += convertion[e->getEndVertex()]->z() * tbf;
    }

    else if ((type==1 && crease_edges.size()>0) || (type==0 && crease_edges.size()>2)){
        //if boundary and with crease edges connected
        x += convertion[e->getEndVertex()]->x();
        y += convertion[e->getEndVertex()]->y();
        z += convertion[e->getEndVertex()]->z();

    }

    return std::make_tuple(x,y,z);





};


// =================================================================
// SIMPLIFICATION
// =================================================================


bool Mesh::isLegalEdge(Edge* e){
    //collect all vertices around the edge
    std::vector<Vertex*> around_vertex = extract_neighbor_vectex(e);
    //insert all vertex into a set to figure if any two of them are the same
    //if any two of them are the same, edge is illegal
    std::set<Vertex*> compare_set;
    for (int i = 0; i < around_vertex.size(); ++i) {
        if(compare_set.insert(around_vertex[i]).second == false){
            return false;
        }
    }
    return true;
}


std::vector<Triangle*> Mesh::extract_triangles(Edge* e){
    //extract all triangles relevant to this edge

    std::vector<Triangle*> result;
    Edge *edge_working_on = e;

    //clockwise search the vertex
    bool all_reached=false;
    if (edge_working_on->getNext()->getOpposite() != NULL) {
        edge_working_on = edge_working_on->getNext()->getOpposite();
        while ((edge_working_on->getNext()->getOpposite() != NULL) &&
               (edge_working_on->getNext()->getOpposite() != e))
        {
            result.push_back(edge_working_on->getTriangle());
            edge_working_on = edge_working_on->getNext()->getOpposite();
        }
        if(edge_working_on->getNext()->getOpposite() != NULL &&
           edge_working_on->getNext()->getOpposite() == e)
        {
            all_reached=true;
        }
        if(edge_working_on->getNext()->getOpposite() == NULL){
            result.push_back(edge_working_on->getTriangle());
        }
    }

    edge_working_on=e;

    //counter clockwise search the vertex if not all edge reached
    if(all_reached==false) {
        if ((edge_working_on->getOpposite() != NULL)&&
            edge_working_on->getOpposite()->getNext()->getNext()->getOpposite()!=NULL)
        {
            edge_working_on = edge_working_on->getOpposite()->getNext()->getNext()->getOpposite()->getNext()->getNext();
            while ((edge_working_on->getOpposite() != NULL)) {
                result.push_back(edge_working_on->getTriangle());
                edge_working_on = edge_working_on->getOpposite()->getNext()->getNext();
            }
            result.push_back(edge_working_on->getTriangle());
        }

    }

    edge_working_on=e;

    //if other side can not be reached by counter clockwise iteration
    if((edge_working_on->getOpposite() == NULL) &&
       (edge_working_on->getNext()->getNext()->getOpposite()!=NULL)){
        edge_working_on = e->getNext()->getNext()->getOpposite()->getNext()->getNext();
        while(edge_working_on->getOpposite()!=NULL){
            result.push_back(edge_working_on->getTriangle());
            edge_working_on = edge_working_on->getOpposite()->getNext()->getNext();
        }
        result.push_back(edge_working_on->getTriangle());
    }


    return result;
}


std::vector<Vertex*> Mesh::extract_half_neighbor_vectex(Edge* e){
    //extract vertices at the end of edge e
    std::vector<Vertex*> result;
    Edge *edge_working_on = e;
    //clockwise search the vertex
    bool all_reached=false;
    if (edge_working_on->getNext()->getOpposite() != NULL) {
        edge_working_on = edge_working_on->getNext()->getOpposite();
        while ((edge_working_on->getNext()->getOpposite() != NULL) &&
               (edge_working_on->getNext()->getOpposite() != e))
        {
            result.push_back(edge_working_on->getStartVertex());
            edge_working_on = edge_working_on->getNext()->getOpposite();
        }
        if(edge_working_on->getNext()->getOpposite() != NULL &&
           edge_working_on->getNext()->getOpposite() == e)
        {
            all_reached=true;
        }
        if(edge_working_on->getNext()->getOpposite() == NULL){
            result.push_back(edge_working_on->getStartVertex());
            result.push_back(edge_working_on->getNext()->getEndVertex());
        }
    }

    edge_working_on=e;
    //counter clockwise search the vertex if not all edge reached
    if(all_reached==false) {
        if ((edge_working_on->getOpposite() != NULL)&&
            edge_working_on->getOpposite()->getNext()->getNext()->getOpposite()!=NULL)
        {
            edge_working_on = edge_working_on->getOpposite()->getNext()->getNext()->getOpposite()->getNext()->getNext();
            while ((edge_working_on->getOpposite() != NULL) &&
                   (edge_working_on!= e)) {
                result.push_back(edge_working_on->getStartVertex());
                edge_working_on = edge_working_on->getOpposite()->getNext()->getNext();
            }
            if(edge_working_on->getOpposite() == NULL){
                result.push_back(edge_working_on->getStartVertex());
            }
        }
    }

    edge_working_on=e;

    //if other side can not be reached by clockwise iteration
    if((edge_working_on->getOpposite() == NULL) &&
            (edge_working_on->getNext()->getNext()->getOpposite()!=NULL)){
        edge_working_on = e->getNext()->getNext()->getOpposite()->getNext()->getNext();
        while(edge_working_on->getOpposite()!=NULL){
            result.push_back(edge_working_on->getStartVertex());
            edge_working_on = edge_working_on->getOpposite()->getNext()->getNext();
        }
        result.push_back(edge_working_on->getStartVertex());
    }


    return result;

}

std::vector<Vertex*> Mesh::extract_neighbor_vectex(Edge* e){
    //extract all vertices around this edge (exclude edge vertex)
    std::vector<Vertex*> results;
    results = extract_half_neighbor_vectex(e);
    if(e->getOpposite()!=NULL){
        std::vector<Vertex*> temp = extract_half_neighbor_vectex(e->getOpposite());
        results.reserve( results.size() + temp.size() );
        results.insert(results.end(),temp.begin(),temp.end());
    }
    return results;
}


void Mesh::Simplification(int target_tri_count) {
    // clear out any previous relationships between vertices
    vertex_parents.clear();
    
    printf ("Simplify the mesh! %d -> %d\n", numTriangles(), target_tri_count);

    while(numTriangles()>target_tri_count) {
        //Find the shortest edge
        float shortest = -1;
        float weight;
        float crease = 0;
        Edge *shortest_edge;
        for (edgeshashtype::iterator iter = edges.begin(); iter != edges.end(); iter++) {
            //if preserve crease is on we calculate weight with crease invloved
            if(args->preserve_crease){
                weight = iter->second->Length();
                weight += iter->second->getNext()->getCrease()/2.0f;
                weight += iter->second->getNext()->getNext()->getCrease()/2.0f;
                weight += iter->second->getCrease();
                if(iter->second->getOpposite()!=NULL){
                    weight += iter->second->getOpposite()->getNext()->getCrease()/2.0f;
                    weight += iter->second->getOpposite()->getNext()->getNext()->getCrease()/2.0f;
                    weight += iter->second->getOpposite()->getCrease()/2.0f;
                }

            }
            else{
                //else use regular short first
                weight = iter->second->Length();

            }

            if (shortest > 0) {
                if (weight < shortest) {
                    if (isLegalEdge(iter->second)) {
                        shortest = weight;
                        shortest_edge = iter->second;
                    }
                }
            } else if (shortest < 0) {
                if (isLegalEdge(iter->second)) {
                    shortest = weight;
                    shortest_edge = iter->second;
                }
            }
        }

        crease = shortest_edge->getCrease();

        //if no legal edge found
        if(shortest==-1){
            std::cout<<"no legal edge found\n";
            return;
        }



        std::vector<Triangle *> to_remove;
        std::vector<Triangle *> to_adjust;
        std::vector<std::tuple<Vertex *, Vertex *, Vertex *>> to_add;
        //Collect triangles to be removed
        to_remove.push_back(shortest_edge->getTriangle());
        if (shortest_edge->getOpposite() != NULL) {
            to_remove.push_back(shortest_edge->getOpposite()->getTriangle());
        }

        //collect triangles
        to_adjust = extract_triangles(shortest_edge);
        if(shortest_edge->getOpposite()!=NULL){
            std::vector<Triangle*> temp = extract_triangles(shortest_edge->getOpposite());
            to_adjust.reserve( to_adjust.size() + temp.size() );
            to_adjust.insert(to_adjust.end(),temp.begin(),temp.end());
        }

        //adjust triangles

        //Vertex operation
        //identify vertex to be removed
        Vertex *vertex_to_remove_1 = shortest_edge->getStartVertex();
        Vertex *vertex_to_remove_2 = shortest_edge->getEndVertex();


        //calculate the new vertex
        double new_x = (vertex_to_remove_1->x() + vertex_to_remove_2->x()) / 2.0f;
        double new_y = (vertex_to_remove_1->y() + vertex_to_remove_2->y()) / 2.0f;
        double new_z = (vertex_to_remove_1->z() + vertex_to_remove_2->z()) / 2.0f;

        //add the new vertex
        Vertex *new_vertex = addVertex(glm::vec3(new_x, new_y, new_z));
        for (int i = 0; i < to_adjust.size(); ++i) {
            //The rest two vertices for reconstruct triangles with new vertex;
            std::tuple<Vertex *, Vertex *, Vertex *> temp;
            //find the vertex not belong to the edge to be remove
            if ((*to_adjust[i])[0] == vertex_to_remove_1 or
                (*to_adjust[i])[0] == vertex_to_remove_2) {
                temp = std::make_tuple(new_vertex, (*to_adjust[i])[1], (*to_adjust[i])[2]);
            } else if ((*to_adjust[i])[1] == vertex_to_remove_1 or
                       (*to_adjust[i])[1] == vertex_to_remove_2) {
                temp = std::make_tuple((*to_adjust[i])[0], new_vertex, (*to_adjust[i])[2]);
            } else if ((*to_adjust[i])[2] == vertex_to_remove_1 or
                    (*to_adjust[i])[2] == vertex_to_remove_2) {
                temp = std::make_tuple((*to_adjust[i])[0], (*to_adjust[i])[1], new_vertex);
            }
            removeTriangle(to_adjust[i]);
            to_add.push_back(temp);
        }


        //remove triangle need to be removed
        assert(to_remove.size()<=2);
        for (int j = 0; j < to_remove.size(); ++j) {
            removeTriangle(to_remove[j]);
        }

        //repair holes on the mesh by adding triangles
        for (int k = 0; k < to_add.size(); ++k) {
            addTriangle(std::get<0>(to_add[k]), std::get<1>(to_add[k]), std::get<2>(to_add[k]));
            if (args->preserve_crease && crease!=0){
                Edge *e1 =  edges[std::make_pair(std::get<0>(to_add[k]), std::get<1>(to_add[k]))];
                Edge *e2 =  edges[std::make_pair(std::get<1>(to_add[k]), std::get<2>(to_add[k]))];
                Edge *e3 =  edges[std::make_pair(std::get<2>(to_add[k]), std::get<0>(to_add[k]))];
                e1->setCrease(crease);
                e2->setCrease(crease);
                e3->setCrease(crease);
                if (e1->getOpposite()!=NULL) {e1->getOpposite()->setCrease(crease);}
                if (e2->getOpposite()!=NULL) {e2->getOpposite()->setCrease(crease);}
                if (e3->getOpposite()!=NULL) {e3->getOpposite()->setCrease(crease);}
            }
        }

    }

}


// =================================================================
