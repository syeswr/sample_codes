#include "glCanvas.h"

#include "raytracer.h"
#include "material.h"
#include "argparser.h"
#include "raytree.h"
#include "utils.h"
#include "mesh.h"
#include "face.h"
#include "primitive.h"
#include "photon_mapping.h"


// ===========================================================================
// casts a single ray through the scene geometry and finds the closest hit
bool RayTracer::CastRay(const Ray &ray, Hit &h, bool use_rasterized_patches) const {
  bool answer = false;

  // intersect each of the quads
  for (int i = 0; i < mesh->numOriginalQuads(); i++) {
    Face *f = mesh->getOriginalQuad(i);
    if (f->intersect(ray,h,args->intersect_backfacing)) answer = true;
  }

  // intersect each of the primitives (either the patches, or the original primitives)
  if (use_rasterized_patches) {
    for (int i = 0; i < mesh->numRasterizedPrimitiveFaces(); i++) {
      Face *f = mesh->getRasterizedPrimitiveFace(i);
      if (f->intersect(ray,h,args->intersect_backfacing)) answer = true;
    }
  } else {
    int num_primitives = mesh->numPrimitives();
    for (int i = 0; i < num_primitives; i++) {
      if (mesh->getPrimitive(i)->intersect(ray,h)) answer = true;
    }
  }
  return answer;
}

// ===========================================================================
// does the recursive (shadow rays & recursive rays) work
glm::vec3 RayTracer::TraceRay(Ray &ray, Hit &hit, int bounce_count) const {

  // First cast a ray and see if we hit anything.
  hit = Hit();
  bool intersect = CastRay(ray,hit,false);
    
  // if there is no intersection, simply return the background color
  if (intersect == false) {
    return glm::vec3(srgb_to_linear(mesh->background_color.r),
                     srgb_to_linear(mesh->background_color.g),
                     srgb_to_linear(mesh->background_color.b));
  }

  // otherwise decide what to do based on the material
  Material *m = hit.getMaterial();
  assert (m != NULL);
  if (glm::length(m->getEmittedColor()) > 0.001) {
    return glm::vec3(1,1,1);
  } 
 
  
  glm::vec3 normal = hit.getNormal();
  glm::vec3 point = ray.pointAtParameter(hit.getT());
  glm::vec3 answer;

  // ----------------------------------------------
  //  start with the indirect light (ambient light)
  glm::vec3 diffuse_color = m->getDiffuseColor(hit.get_s(),hit.get_t());
  if (args->gather_indirect) {
    // photon mapping for more accurate indirect light
    answer = diffuse_color * (photon_mapping->GatherIndirect(point, normal, ray.getDirection()) + args->ambient_light);
  } else {
    // the usual ray tracing hack for indirect light
    answer = diffuse_color * args->ambient_light;
  }      

  // ----------------------------------------------
  // add contributions from each light that is not in shadow
  int num_lights = mesh->getLights().size();
  for (int i = 0; i < num_lights; i++) {

    Face *f = mesh->getLights()[i];
    glm::vec3 lightColor = f->getMaterial()->getEmittedColor() * f->getArea();
    glm::vec3 myLightColor;
    glm::vec3 lightCentroid = f->computeCentroid();
    glm::vec3 dirToLightCentroid = glm::normalize(lightCentroid-point);
    

    // ===========================================
    I AM GOING TO ADD SHADOW & SOFT SHADOW LOGIC HERE
    // ===========================================

    float distToLightCentroid = glm::length(lightCentroid-point);
    myLightColor = lightColor / float(M_PI*distToLightCentroid*distToLightCentroid);




    
    // add the lighting contribution from this particular light at this point
    answer += m->Shade(ray,hit,dirToLightCentroid,myLightColor,args);
  }
      
  // ----------------------------------------------
  // add contribution from reflection, if the surface is shiny
  glm::vec3 reflectiveColor = m->getReflectiveColor();

  if (bounce_count > 0) {
    glm::vec3 n = glm::normalize(normal);
    Ray r_ray = Ray(point, ray.getDirection() - 2.0f * glm::dot(ray.getDirection(), n) * n);
    answer += reflectiveColor * TraceRay(r_ray, hit, bounce_count - 1);
  }

  
  return answer; 

}







bool Sphere::intersect(const Ray &r, Hit &h) const {



    float answer_t;
#pragma clang diagnostic push
#pragma ide diagnostic ignored "IncompatibleTypes"
    float a = glm::dot(r.getDirection(), r.getDirection());
#pragma clang diagnostic pop
#pragma clang diagnostic push
#pragma ide diagnostic ignored "IncompatibleTypes"
    float b = glm::dot(r.getDirection(), 2.0f * r.getOrigin()) -
            glm::dot( 2.0f * r.getDirection(), center );
#pragma clang diagnostic pop

#pragma clang diagnostic push
#pragma ide diagnostic ignored "IncompatibleTypes"
  float c = glm::dot(r.getOrigin(), r.getOrigin()) -
                2.0f * glm::dot( r.getOrigin(), center) +
                glm::dot(center, center) -
                radius * radius;
#pragma clang diagnostic pop
    glm::vec3 normal;
    float delta = b*b - 4.0f * a * c;
    if  ( delta <= 0 )
    {
        return false;
    }
    else{
        //we want the smaller t, which is the first intersection
        //between the light and the sphere
        float answer_t_1;
        float answer_t_2;
        float d = sqrt(delta);
        answer_t_1 = (-b + d) / (2.0f * a);
        answer_t_2 = (-b - d) / (2.0f * a);
        if (answer_t_1 > 0 && answer_t_2 > 0) {
            answer_t = answer_t_1 > answer_t_2 ? answer_t_2 : answer_t_1;
        }
        else if (answer_t_1 > 0.001 && answer_t_2 <= 0){
            answer_t = answer_t_1;
        }
        else if (answer_t_2 > 0.001 && answer_t_1 <= 0){
            answer_t = answer_t_2;
        }
        else{
            return false;
        }

        glm::vec3 intersection_p = r.getOrigin() + r.getDirection() * answer_t;
        normal = glm::normalize(intersection_p - center);
        if (glm::distance(r.getOrigin(), center) < radius){
            normal = -1.0f * normal;
        }
    }

    h.set(answer_t, this->getMaterial(), normal);


  // return true if the sphere was intersected, and update the hit
  // data structure to contain the value of t for the ray at the
  // intersection point, the material, and the normal


  return true;

} 
