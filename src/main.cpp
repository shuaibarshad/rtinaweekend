#include <iostream>
#include "main.h"
#include "sphere.h"
#include "moving_sphere.h"
#include "hitable_list.h"
#include "float.h"
#include "camera.h"
#include "material.h"
#include "bvh_node.h"
#include "texture.h"
#include "xy_rect.h"

#if defined(WIN32) or defined(LINUX)
#define MAXFLOAT FLT_MAX
#endif

#define NUM_BOUNCE 5
#define SAMPLES_PER_PIXEL 16
//#define USE_BVH

unsigned long long ray_prim_intersect_count = 0;
unsigned long long ray_node_intersect_count = 0; 


vec3 color(const ray& r, hitable* world, int depth) {
    hit_record rec;
    if (world->hit(r, 0.001, MAXFLOAT, rec)) {
        ray scattered;
        vec3 attenuation;
        vec3 emitted = rec.mat_ptr->emitted(rec.u, rec.v, rec.p);
        if (depth < NUM_BOUNCE && rec.mat_ptr->scatter(r, rec, attenuation, scattered)) {
            return attenuation*color(scattered, world, depth+1);
        }
        else {
            return vec3(0, 0, 0);
        }
    }
    // Commenting out to test rectangle lights
    return vec3(0, 0, 0);
    //else {
    //    vec3 unit_direction = unit_vector(r.direction());
    //    float t = 0.5*(unit_direction.y() + 1.0);
    //    return (1.0-t)*vec3(1.0, 1.0, 1.0) + t*vec3(0.5, 0.7, 1.0);
    //}
}

hitable* random_scene() {
    int n = 500;
    hitable** list = new hitable*[n+1];
    texture* checker = new checker_texture(
            new constant_texture(vec3(0.2, 0.2, 0.1)), 
            new constant_texture(vec3(0.9, 0.9, 0.9))
            );
    list[0] = new sphere(vec3(0, -1000, 0), 1000, new lambertian(checker));
    //list[0] = new sphere(vec3(0, -1000, 0), 1000, new lambertian(new constant_texture(vec3(0.5, 0.5, 0.5))));
    int i = 1;
    for (int a = -11; a < 11; a++) {
        for (int b = -11; b < 11; b++) {
            float choose_mat = drand48();
            vec3 center(a+0.9*drand48(), 0.2, b+0.9*drand48());
            if ((center-vec3(4,0.2,0)).length() > 0.9) {
                if (choose_mat < 0.8) {  // diffuse
                    list[i++] = new sphere(center, 0.2, 
                            new lambertian(new constant_texture(vec3(drand48()*drand48(), drand48()*drand48(), drand48()*drand48()))));
                }
                else if(choose_mat < 0.95) {  // metal
                    list[i++] = new sphere(center, 0.2, 
                            new metal(new constant_texture(vec3(0.5*(1+drand48()), 0.5*(1+drand48()), 0.5*(1+drand48()))), 0.5*drand48()));
                }
                else {  // glass
                    list[i++] = new sphere(center, 0.2, new dielectric(1.5));
                }
            }
        }
    }

    list[i++] = new sphere(vec3(0,1,0), 1.0, new dielectric(1.5));
    list[i++] = new sphere(vec3(-4,1,0), 1.0, new lambertian(new constant_texture(vec3(0.4,0.2,0.1))));
    list[i++] = new sphere(vec3(4,1,0), 1.0, new metal(new constant_texture(vec3(0.7, 0.6, 0.5)), 0.0));

    return new hitable_list(list, i);
}

hitable* random_scene_moving_spheres() {
    int n = 50000;
    hitable** list = new hitable*[n+1];
    list[0] = new sphere(vec3(0, -1000, 0), 1000, new lambertian(new constant_texture(vec3(0.5, 0.5, 0.5))));
    int i = 1;
    for (int a = -10; a < 10; a++) {
        for (int b = -10; b < 10; b++) {
            float choose_mat = drand48();
            vec3 center(a+0.9*drand48(), 0.2, b+0.9*drand48());
            if ((center-vec3(4,0.2,0)).length() > 0.9) {
                if (choose_mat < 0.8) {       // diffuse
                    list[i++] = new moving_sphere(center, center+vec3(0, 0.5*drand48(), 0), 0.0, 1.0, 0.2, 
                            new lambertian(new constant_texture(vec3(drand48()*drand48(), drand48()*drand48(), drand48()*drand48()))));
                }
                else if(choose_mat < 0.95) {  // metal
                    list[i++] = new sphere(center, 0.2, 
                            new metal(new constant_texture(vec3(0.5*(1+drand48()), 0.5*(1+drand48()), 0.5*(1+drand48()))), 0.5*drand48()));
                }
                else {  // glass
                    list[i++] = new sphere(center, 0.2, new dielectric(1.5));
                }
            }
        }
    }

    list[i++] = new sphere(vec3(0,1,0), 1.0, new dielectric(1.5));
    list[i++] = new sphere(vec3(-4,1,0), 1.0, new lambertian(new constant_texture(vec3(0.4,0.2,0.1))));
    list[i++] = new sphere(vec3(4,1,0), 1.0, new metal(new constant_texture(vec3(0.7, 0.6, 0.5)), 0.0));

    return new hitable_list(list, i);
}

hitable* two_perlin_spheres() {
    texture* pertext = new noise_texture(1);
    hitable** list = new hitable*[2];
    list[0] = new sphere(vec3(0,-1000,0), 1000, new lambertian(pertext));
    list[1] = new sphere(vec3(0,2,0), 2, new lambertian(pertext));
    return new hitable_list(list, 2);
}

hitable* earth() {
    int nx, ny, nn;
    unsigned char* tex_data = stbi_load("earthmap.jpg", &nx, &ny, &nn, 0);
	if (!tex_data) {
		std::cerr << "Couldn't load image texture";
		exit(-1);
	}
    material* mat = new lambertian((texture*)new image_texture(tex_data, nx, ny));
    hitable** list = new hitable*[1];
    list[0] = new sphere(vec3(0,0,0), 3, mat);
    return new hitable_list(list, 1);
}

hitable* simple_light() {
    texture* pertext = new noise_texture(4);
    hitable** list = new hitable*[4];
    list[0] = new sphere(vec3(0, -1000, 0), 1000, new lambertian(pertext));
    list[1] = new sphere(vec3(0, 2, 0), 2, new lambertian(pertext));
    list[2] = new sphere(vec3(0, 7, 0), 2, new diffuse_light(new constant_texture(vec3(4,4,4))));
    list[3] = new xy_rect(3, 5, 1, 3, -2, new diffuse_light(new constant_texture(vec3(4,4,4))));
    return new hitable_list(list, 4);
}

int main(int argc, char* argv[])
{
    int nx = 800;
    int ny = 400;
    int ns = SAMPLES_PER_PIXEL;

    std::cout << "P3\n" << nx << " " << ny << "\n255\n";

#if 0
    hitable* list[4];
    list[0] = new sphere(vec3(0,0,-1), 0.5, new lambertian(new constant_texture(vec3(0.1, 0.2, 0.5))));
    list[1] = new sphere(vec3(0,-100.5,-1), 100, new lambertian(new constant_texture(vec3(0.8, 0.8, 0.0))));
    list[2] = new sphere(vec3(1,0,-1), 0.5, new metal(vec3(0.8, 0.6, 0.2)));
    list[3] = new sphere(vec3(-1,0,-1), 0.5, new dielectric(1.5));
    hitable* world = new hitable_list(list, 4);
    vec3 lookfrom(3,3,2);
    vec3 lookat(0,0,-1);
    float dist_to_focus = (lookfrom - lookat).length();
    float aperture = 2.0;
    camera cam(lookfrom, lookat, vec3(0,1,0), 45, float(nx)/float(ny), aperture, dist_to_focus, 0.0, 1.0);
#endif
    
    float t0 = 0.0;
    float t1 = 1.0;
    //hitable* world = random_scene();
    //hitable* world = random_scene_moving_spheres();
    //vec3 lookfrom(8,1.5,4);
    //// Geometry and camera for Perlin noise
    //hitable* world = two_perlin_spheres();
    // Earth 
    hitable* world = earth();

#ifdef USE_BVH
    bvh_node* world_bvh = new bvh_node(((hitable_list*)world)->list, ((hitable_list*)world)->list_size, t0,  t1);
#endif
    vec3 lookfrom(13,2,3);
    vec3 lookat(0,0,0);
    float dist_to_focus = 10;
    //float aperture = 0.1;
    float aperture = 0.0;
    camera cam(lookfrom, lookat, vec3(0,1,0), 20, float(nx)/float(ny), aperture, dist_to_focus, t0, t1);
    //// Camera location for the earth sphere
    //camera cam(lookfrom, lookat, vec3(0,1,0), 45, float(nx)/float(ny), aperture, dist_to_focus, t0, t1);

    //float R = cos(M_PI/4);
    //hitable* list[2];
    //list[0] = new sphere(vec3(-R, 0, -1), R, new lambertian(new constant_texture(vec3(0, 0, 1))));
    //list[1] = new sphere(vec3( R, 0, -1), R, new lambertian(new constant_texture(vec3(1, 0, 0))));
    //hitable* world = new hitable_list(list, 2);
    //camera cam(90, float(nx)/float(ny));

    for (int j = ny-1; j >= 0; j--) {
        for (int i = 0; i < nx; i++) {
            vec3 col(0, 0, 0);
            for (int s = 0; s < ns; s++) {
                float u = float(i + drand48()) / float(nx);
                float v = float(j + drand48()) / float(ny);
                ray r = cam.get_ray(u, v);
                vec3 p = r.point_at_parameter(2.0);
#ifdef USE_BVH
				col += color(r, world_bvh, 0);
#else
                col += color(r, world, 0);
#endif
            }
            col /= float(ns);
            col = vec3(sqrt(col[0]), sqrt(col[1]), sqrt(col[2]));
            int ir = int(255.99*col[0]);
            int ig = int(255.99*col[1]);
            int ib = int(255.99*col[2]);

            std::cout << ir << " " << ig << " " << ib << "\n";
        }
    }

    std::cerr << "ray-prim tests: " << ray_prim_intersect_count << "\n";
    std::cerr << "ray-node tests: " << ray_node_intersect_count << "\n";
}

