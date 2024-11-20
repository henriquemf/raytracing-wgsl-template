const THREAD_COUNT = 16;
const RAY_TMIN = 0.0001;
const RAY_TMAX = 100.0;
const PI = 3.1415927f;
const FRAC_1_PI = 0.31830987f;
const FRAC_2_PI = 1.5707964f;

@group(0) @binding(0)  
  var<storage, read_write> fb : array<vec4f>;

@group(0) @binding(1)
  var<storage, read_write> rtfb : array<vec4f>;

@group(1) @binding(0)
  var<storage, read_write> uniforms : array<f32>;

@group(2) @binding(0)
  var<storage, read_write> spheresb : array<sphere>;

@group(2) @binding(1)
  var<storage, read_write> quadsb : array<quad>;

@group(2) @binding(2)
  var<storage, read_write> boxesb : array<box>;

@group(2) @binding(3)
  var<storage, read_write> trianglesb : array<triangle>;

@group(2) @binding(4)
  var<storage, read_write> meshb : array<mesh>;

struct ray {
  origin : vec3f,
  direction : vec3f,
};

struct sphere {
  transform : vec4f,
  color : vec4f,
  material : vec4f,
};

struct quad {
  Q : vec4f,
  u : vec4f,
  v : vec4f,
  color : vec4f,
  material : vec4f,
};

struct box {
  center : vec4f,
  radius : vec4f,
  rotation: vec4f,
  color : vec4f,
  material : vec4f,
};

struct triangle {
  v0 : vec4f,
  v1 : vec4f,
  v2 : vec4f,
};

struct mesh {
  transform : vec4f,
  scale : vec4f,
  rotation : vec4f,
  color : vec4f,
  material : vec4f,
  min : vec4f,
  max : vec4f,
  show_bb : f32,
  start : f32,
  end : f32,
};

struct material_behaviour {
  scatter : bool,
  direction : vec3f,
};

struct camera {
  origin : vec3f,
  lower_left_corner : vec3f,
  horizontal : vec3f,
  vertical : vec3f,
  u : vec3f,
  v : vec3f,
  w : vec3f,
  lens_radius : f32,
};

struct hit_record {
  t : f32,
  p : vec3f,
  normal : vec3f,
  object_color : vec4f,
  object_material : vec4f,
  frontface : bool,
  hit_anything : bool,
};

fn ray_at(r: ray, t: f32) -> vec3f
{
  return r.origin + t * r.direction;
}

fn get_ray(cam: camera, uv: vec2f, rng_state: ptr<function, u32>) -> ray
{
  var rd = cam.lens_radius * rng_next_vec3_in_unit_disk(rng_state);
  var offset = cam.u * rd.x + cam.v * rd.y;
  return ray(cam.origin + offset, normalize(cam.lower_left_corner + uv.x * cam.horizontal + uv.y * cam.vertical - cam.origin - offset));
}

fn get_camera(lookfrom: vec3f, lookat: vec3f, vup: vec3f, vfov: f32, aspect_ratio: f32, aperture: f32, focus_dist: f32) -> camera
{
  var camera = camera();
  camera.lens_radius = aperture / 2.0;

  var theta = degrees_to_radians(vfov);
  var h = tan(theta / 2.0);
  var w = aspect_ratio * h;

  camera.origin = lookfrom;
  camera.w = normalize(lookfrom - lookat);
  camera.u = normalize(cross(vup, camera.w));
  camera.v = cross(camera.u, camera.w);

  camera.lower_left_corner = camera.origin - w * focus_dist * camera.u - h * focus_dist * camera.v - focus_dist * camera.w;
  camera.horizontal = 2.0 * w * focus_dist * camera.u;
  camera.vertical = 2.0 * h * focus_dist * camera.v;

  return camera;
}

fn environment_color(direction: vec3f, color1: vec3f, color2: vec3f) -> vec3f
{
  var unit_direction = normalize(direction);
  var t = 0.5 * (unit_direction.y + 1.0);
  var col = (1.0 - t) * color1 + t * color2;

  var sun_direction = normalize(vec3(uniforms[13], uniforms[14], uniforms[15]));
  var sun_color = int_to_rgb(i32(uniforms[17]));
  var sun_intensity = uniforms[16];
  var sun_size = uniforms[18];

  var sun = clamp(dot(sun_direction, unit_direction), 0.0, 1.0);
  col += sun_color * max(0, (pow(sun, sun_size) * sun_intensity));

  return col;
}

fn check_ray_collision(r: ray, max: f32) -> hit_record
{
  var spheresCount = i32(uniforms[19]);
  var quadsCount = i32(uniforms[20]);
  var boxesCount = i32(uniforms[21]);
  var trianglesCount = i32(uniforms[22]);
  var meshCount = i32(uniforms[27]);

  var record = hit_record(RAY_TMAX, vec3f(0.0), vec3f(0.0), vec4f(0.0), vec4f(0.0), false, false);
  var closest = record;

  for (var i = 0; i < spheresCount; i = i + 1) {
    var sphere = spheresb[i];
    hit_sphere(sphere.transform.xyz, sphere.transform.w, r, &record, max);
    
    if (record.hit_anything && record.t < closest.t) {
        closest = record;
        closest.object_color = sphere.color;
        closest.object_material = sphere.material;
    }
  }

  // Check for sphere collisions
    for (var i = 0; i < spheresCount; i = i + 1) {
        var sphere = spheresb[i];
        var temp_record = hit_record(RAY_TMAX, vec3f(0.0), vec3f(0.0), vec4f(0.0), vec4f(0.0), false, false);
        hit_sphere(sphere.transform.xyz, sphere.transform.w, r, &temp_record, max);
        
        if (temp_record.hit_anything && temp_record.t < closest.t) {
            closest = temp_record;
            closest.object_color = sphere.color;
            closest.object_material = sphere.material;
        }
    }

    // Check for quad collisions
    for (var i = 0; i < quadsCount; i = i + 1) {
        var quad = quadsb[i];
        var temp_record = hit_record(RAY_TMAX, vec3f(0.0), vec3f(0.0), quad.color, quad.material, false, true);
        hit_quad(r, quad.Q, quad.u, quad.v, &temp_record, max);
        
        if (temp_record.hit_anything && temp_record.t < closest.t) {
            closest = temp_record;
            closest.object_color = quad.color;
            closest.object_material = quad.material;
        }
    }

    // Check for box collisions
    for (var i = 0; i < boxesCount; i = i + 1) {
        var box = boxesb[i];
        var temp_record = hit_record(RAY_TMAX, vec3f(0.0), vec3f(0.0), box.color, box.material, false, true);
        hit_box(r, box.center.xyz, box.radius.xyz, &temp_record, max);
        
        if (temp_record.hit_anything && temp_record.t < closest.t) {
            closest = temp_record;
            closest.object_color = box.color;
            closest.object_material = box.material;
        }
    }

    // Check for triangle collisions
    for (var i = 0; i < trianglesCount; i = i + 1) {
        var triangle = trianglesb[i];
        var temp_record = hit_record(RAY_TMAX, vec3f(0.0), vec3f(0.0), vec4f(1.0), vec4f(1.0), false, true);
        hit_triangle(r, triangle.v0.xyz, triangle.v1.xyz, triangle.v2.xyz, &temp_record, max);
        
        if (temp_record.hit_anything && temp_record.t < closest.t) {
            closest = temp_record;
            closest.object_color = vec4f(1.0); // Supondo uma cor padrão para triângulos
            closest.object_material = vec4f(1.0); // Supondo um material padrão
        }
    }

  return closest;
}

fn lambertian(normal : vec3f, absorption: f32, random_sphere: vec3f, rng_state: ptr<function, u32>) -> material_behaviour
{
  var prob = rng_next_float(rng_state);
  if (prob < absorption) {
      return material_behaviour(false, vec3f(0.0));
  } else {
      var scatter_direction = normal + random_sphere;
      return material_behaviour(true, scatter_direction);
  }
}

fn metal(normal : vec3f, direction: vec3f, fuzz: f32, random_sphere: vec3f) -> material_behaviour
{
  var reflected = direction - 2.0 * dot(direction, normal) * normal;
  var scatter_direction = reflected + fuzz * random_sphere;
  return material_behaviour(true, normalize(scatter_direction));
}

fn dielectric(normal : vec3f, r_direction: vec3f, refraction_index: f32, frontface: bool, random_sphere: vec3f, fuzz: f32, rng_state: ptr<function, u32>) -> material_behaviour
{  
  var refraction_ratio = select(1.0 / refraction_index, refraction_index, frontface); // Ratio of indices of refraction
  var cos_theta = min(dot(-r_direction, normal), 1.0); // Angle of incidence
  var sin_theta = sqrt(1.0 - cos_theta * cos_theta);

  // Check for total internal reflection
  if (refraction_ratio * sin_theta > 1.0) {
      // Reflect instead of refracting
      var reflected = r_direction - 2.0 * dot(r_direction, normal) * normal;
      return material_behaviour(true, normalize(reflected));
  }

  // Calculate Fresnel reflectance (Schlick's approximation)
  var r0 = pow((1.0 - refraction_ratio) / (1.0 + refraction_ratio), 2.0);
  var reflectance = r0 + (1.0 - r0) * pow(1.0 - cos_theta, 5.0);

  // Randomly reflect or refract based on reflectance
  if (rng_next_float(rng_state) < reflectance) {
      var reflected = r_direction - 2.0 * dot(r_direction, normal) * normal;
      return material_behaviour(true, normalize(reflected));
  } else {
      // Refract the ray
      var refracted_perpendicular = refraction_ratio * (r_direction + cos_theta * normal);
      var refracted_parallel = -sqrt(abs(1.0 - dot(refracted_perpendicular, refracted_perpendicular))) * normal;
      var refracted = refracted_perpendicular + refracted_parallel;
      return material_behaviour(true, refracted);
  }
}

fn emmisive(color: vec3f, light: f32) -> material_behaviour
{
    return material_behaviour(false, color * light);
}

fn trace(r: ray, rng_state: ptr<function, u32>) -> vec3f {
    var max_bounces = i32(uniforms[2]); // Máximo de reflexões
    var accumulated_light = vec3f(0.0); // Luz acumulada
    var current_color = vec3f(1.0); // Multiplicador de cor atual
    var current_ray = r; // Raio atual

    var bg_color1 = int_to_rgb(i32(uniforms[11])); // Cor de fundo 1
    var bg_color2 = int_to_rgb(i32(uniforms[12])); // Cor de fundo 2
    var material_behaviour = material_behaviour(true, vec3f(0.0)); // Comportamento do material

    for (var j = 0; j < max_bounces; j = j + 1) {
        var hit = check_ray_collision(current_ray, RAY_TMAX);        
        if (!hit.hit_anything) {
            accumulated_light += current_color * environment_color(current_ray.direction, bg_color1, bg_color2);
            break;
        }

        var smoothness = hit.object_material.x;
        var absorption = hit.object_material.y;
        var specular = hit.object_material.z;
        var emission = hit.object_material.w;

        if (emission > 0.0) {
          var emit_behaviour = emmisive(hit.object_color.rgb, emission);
          accumulated_light += current_color * emit_behaviour.direction;
        }
        
        var random_in_sphere = rng_next_vec3_in_unit_sphere(rng_state);
        var random_value = rng_next_float(rng_state);

        if (smoothness > 0.0) {
            // Metal
            if (specular > random_value) {
                material_behaviour = metal(hit.normal, current_ray.direction, absorption, random_in_sphere);
            } else {
                material_behaviour = lambertian(hit.normal, absorption, random_in_sphere, rng_state);
                current_color *= hit.object_color.rgb * (1.0 - absorption);
            }
        } else if (smoothness < 0.0) {
            // Dielectric
            material_behaviour = dielectric(hit.normal, current_ray.direction, specular, hit.frontface, random_in_sphere, absorption, rng_state);
            current_ray = ray(hit.p, material_behaviour.direction);
            continue;
        } else {
            // Lambertian
            material_behaviour = lambertian(hit.normal, absorption, random_in_sphere, rng_state);
            current_color *= hit.object_color.rgb * (1.0 - absorption);
        }


        if (!material_behaviour.scatter) {
            break;
        }

        let offset = hit.normal * 0.0001;
        current_ray = ray(hit.p + offset, material_behaviour.direction);
    }

    return accumulated_light;
}

@compute @workgroup_size(THREAD_COUNT, THREAD_COUNT, 1)
fn render(@builtin(global_invocation_id) id : vec3u)
{
    var rez = uniforms[1];
    var time = u32(uniforms[0]);

    // init_rng (random number generator) we pass the pixel position, resolution and frame
    var rng_state = init_rng(vec2(id.x, id.y), vec2(u32(rez)), time);

    // Get uv
    var fragCoord = vec2f(f32(id.x), f32(id.y));
    var uv = (fragCoord + sample_square(&rng_state)) / vec2(rez);

    // Camera setup
    var lookfrom = vec3(uniforms[7], uniforms[8], uniforms[9]);
    var lookat = vec3(uniforms[23], uniforms[24], uniforms[25]);
    var cam = get_camera(lookfrom, lookat, vec3(0.0, 1.0, 0.0), uniforms[10], 1.0, uniforms[6], uniforms[5]);
    
    // Get number of samples per pixel from uniforms
    var samples_per_pixel = i32(uniforms[4]);
    
    // Initialize accumulated color
    var color = vec3f(0.0);
    
    // 1. Loop for each sample per pixel
    for(var s = 0; s < samples_per_pixel; s++) {
        // 2. Get ray for this sample
        var r = get_ray(cam, uv, &rng_state);
        
        // 3. Call trace function to get color for this ray
        var sample_color = trace(r, &rng_state);
        
        // Accumulate the color
        color += sample_color;
    }
    
    // 4. Average the color by dividing by number of samples
    color = color / f32(samples_per_pixel);

    var color_out = vec4(linear_to_gamma(color), 1.0);
    
    var map_fb = mapfb(id.xy, rez);
    
    var should_accumulate = uniforms[3];
    var accumulated_color = rtfb[map_fb] * should_accumulate + color_out;

    rtfb[map_fb] = accumulated_color;
    fb[map_fb] = accumulated_color/accumulated_color.w;
}