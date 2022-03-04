use crate::math::Vector3D;
use pixels::{Error, Pixels, SurfaceTexture};
use winit::dpi::LogicalSize;
use winit::event::{Event, VirtualKeyCode};
use winit::event_loop::{ControlFlow, EventLoop};
use winit::window::{WindowBuilder, Window};
use winit_input_helper::WinitInputHelper;
use log::error;
use std::fs::File;
use std::io::BufReader;
use obj::{Obj};

pub const WIDTH: u32 = 200;
pub const HEIGHT: u32 = 200;
pub const BOX_SIZE: i16 = 64;

// #[test]
// fn draw_obj_file_test() {
//     let mut obj: Obj = Obj::load("./obj/african_head.obj").unwrap();
//     let mut obj_buffer = Vec::new();
//     obj.data.write_to_buf(&mut obj_buffer).unwrap();
//     let obj_data = obj::ObjData::load_buf(obj_buffer.as_slice()).unwrap();
//     assert_eq!(obj_data, obj.data);
// }

struct World {
    box_x: i16,
    box_y: i16,
    velocity_x: i16,
    velocity_y: i16,
    should_debug: bool
}

mod rasterizer;

mod math;

fn main() -> Result<(), Error> {
    env_logger::init();
    let event_loop = EventLoop::new();
    let mut input = WinitInputHelper::new();
    let window: Window = {
        let size = LogicalSize::new(WIDTH as f64, HEIGHT as f64);
        WindowBuilder::new()
            .with_title("BFLabs Renderer")
            .with_inner_size(size)
            .with_min_inner_size(size)
            .build(&event_loop)
            .unwrap()
    };

    let mut pixels: Pixels = {
        let window_size = window.inner_size();
        let surface_texture = SurfaceTexture::new(window_size.width,
            window_size.height, &window);
        Pixels::new(WIDTH, HEIGHT, surface_texture).unwrap()
    };
    
    let mut world: World = World::new();
    
    let mut obj: Obj = Obj::load("./obj/african_head.obj").unwrap();

    event_loop.run(move |event, _, control_flow| {
        if let Event::RedrawRequested(_) = event {
            let triangle_vertices: [Vector3D; 3] = [Vector3D::new(10., 10., 0.),
                Vector3D::new(100., 30., 0.), Vector3D::new(190., 160., 0.)];
            // world.draw_triangle(pixels.get_frame(), &triangle_vertices);
            world.draw_obj_file(pixels.get_frame(), &obj);
            if pixels
                .render()
                .map_err(|e| error!("pixels.render() failed: {:?}", e))
                .is_err()
            {
                *control_flow = ControlFlow::Exit;
                return;
            }
        };
        
        if input.update(&event) {
            if input.key_pressed(VirtualKeyCode::Escape) || input.quit() {
                *control_flow = ControlFlow::Exit;
            }

            if let Some(size) = input.window_resized() {
                pixels.resize_surface(size.width, size.height);
            }

            // world.update();
            window.request_redraw();
        };
    });
}

impl World {
    /// Create a new `World` instance that can draw a moving box.
    fn new() -> Self {
        Self {
            box_x: 24,
            box_y: 16,
            velocity_x: 1,
            velocity_y: 1,
            should_debug: true
        }
    }

    /// Update the `World` internal state; bounce the box around the screen.
    fn update(&mut self) {
        if self.box_x <= 0 || self.box_x + BOX_SIZE > WIDTH as i16 {
            self.velocity_x *= -1;
        }
        if self.box_y <= 0 || self.box_y + BOX_SIZE > HEIGHT as i16 {
            self.velocity_y *= -1;
        }

        self.box_x += self.velocity_x;
        self.box_y += self.velocity_y;
    }

    /// Draw the `World` state to the frame buffer.
    ///
    /// Assumes the default texture format: `wgpu::TextureFormat::Rgba8UnormSrgb`
    fn draw_triangle(&self, frame: &mut [u8], pts: &[math::Vector3D; 3]) {
        for (i, pixel) in frame.chunks_exact_mut(4).enumerate() {
            let mut bbox_min: math::Vector3D = math::Vector3D::new((WIDTH - 1) as f32,
                (HEIGHT - 1) as f32, 0.);
            let mut bbox_max: math::Vector3D = math::Vector3D::new(0., 0., 0.);
            let clamp: Vector3D = Vector3D::new(WIDTH as f32 - 1.0,
                HEIGHT as f32 - 1.0, 0.);

            for i in 0..2 {
                bbox_min.x = (0.0 as f32).max(bbox_min.x.min(pts[i].x));
                bbox_min.y = (0.0 as f32).max(bbox_min.y.min(pts[i].y));

                bbox_max.x = clamp.x.min(bbox_max.x.max(pts[i].x));
                bbox_max.y = clamp.y.min(bbox_max.y.max(pts[i].y));
            }

            let x = (i % WIDTH as usize) as i16;
            let y = (WIDTH as usize - i / WIDTH as usize) as i16;

            let temp_point = Vector3D::new(x as f32,
                y as f32, 0.0);
            let bc_screen = Vector3D::barycentric(pts, &temp_point);

            let outside_the_triangle: bool = bc_screen.x < 0.0
                || bc_screen.y < 0.0 || bc_screen.z < 0.0;

            let rgba = if outside_the_triangle {
                [0x5e, 0x48, 0xe8, 0xff]
            } else {
                [0x48, 0xb2, 0xe8, 0xff]
            };

            pixel.copy_from_slice(&rgba);
        }
    }

    fn draw_obj_file(&mut self, frame: &mut [u8], obj: &Obj) {

        let mut should_debug: bool = true;
        
        let mut triple_vertex_iterator: usize = 0;

        // unil the third-to-last element in the vertex array
        while triple_vertex_iterator != obj.data.position.len() - 1 - 3 {
            
            if (self.should_debug)
            {
                println!("{:?}", triple_vertex_iterator);
            }
            let mut triangle_vertices: [Vector3D; 3] = [Vector3D::default(),
                Vector3D::default(), Vector3D::default()];
            // iterate over three vertices 
            let mut triangle_vertices_current_index: usize = 0;
            for i in triple_vertex_iterator..triple_vertex_iterator + 2 {
                let mut single_vertex: Vector3D = Vector3D::default();

                // iterate over items of a vertex
                for j in 0..2 {
                    single_vertex[j] = obj.data.position[i][j];

                    match j {
                        // for X coordinates
                        0 => single_vertex[j] = single_vertex[j]
                            * WIDTH as f32 / 2.,
                        // for Y coordinates
                        1 => single_vertex[j] = single_vertex[j]
                            * HEIGHT as f32 / 2.,
                        // do nothing
                        2 => (),
                        _ => panic!()
                    }
                }

                triangle_vertices[triangle_vertices_current_index] =
                    single_vertex;
                triangle_vertices_current_index =
                    triangle_vertices_current_index + 1;
            }

            self.draw_triangle(frame, &triangle_vertices);
            // let temp_b: bool = triple_vertex_iterator + 3 == obj.data.position.len() - 4;
            // if (temp_b)
            // {
            //     println!("{:?}", "reached end");
            //     should_debug = false;
            //     break;
            // }
            triple_vertex_iterator += 3;
        }
        
        self.should_debug = false;
    }
        
        // for (i, pixel) in frame.chunks_exact_mut(4).enumerate() {
        //     let mut bbox_min: math::Vector3D = math::Vector3D::new((WIDTH - 1) as f32,
        //         (HEIGHT - 1) as f32, 0.);
        //     let mut bbox_max: math::Vector3D = math::Vector3D::new(0., 0., 0.);
        //     let clamp: Vector3D = Vector3D::new(WIDTH as f32 - 1.0,
        //         HEIGHT as f32 - 1.0, 0.);
        // 
        //     for i in 0..2 {
        //         bbox_min.x = (0.0 as f32).max(bbox_min.x.min(pts[i].x));
        //         bbox_min.y = (0.0 as f32).max(bbox_min.y.min(pts[i].y));
        // 
        //         bbox_max.x = clamp.x.min(bbox_max.x.max(pts[i].x));
        //         bbox_max.y = clamp.y.min(bbox_max.y.max(pts[i].y));
        //     }
        // 
        //     let x = (i % WIDTH as usize) as i16;
        //     let y = (WIDTH as usize - i / WIDTH as usize) as i16;
        // 
        //     let temp_point = Vector3D::new(x as f32,
        //         y as f32, 0.0);
        //     let bc_screen = Vector3D::barycentric(pts, &temp_point);
        // 
        //     let outside_the_triangle: bool = bc_screen.x < 0.0
        //         || bc_screen.y < 0.0 || bc_screen.z < 0.0;
        // 
        //     let rgba = if outside_the_triangle {
        //         [0x5e, 0x48, 0xe8, 0xff]
        //     } else {
        //         [0x48, 0xb2, 0xe8, 0xff]
        //     };
        // 
        //     pixel.copy_from_slice(&rgba);
        // }
}