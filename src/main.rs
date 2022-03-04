use crate::math::Vector3D;
use pixels::{Error, Pixels, SurfaceTexture};
use winit::dpi::LogicalSize;
use winit::event::{Event, VirtualKeyCode};
use winit::event_loop::{ControlFlow, EventLoop};
use winit::window::{WindowBuilder, Window};
use winit_input_helper::WinitInputHelper;
use log::error;

pub const WIDTH: u32 = 200;
pub const HEIGHT: u32 = 200;
pub const BOX_SIZE: i16 = 64;

struct World {
    box_x: i16,
    box_y: i16,
    velocity_x: i16,
    velocity_y: i16
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
    
    event_loop.run(move |event, _, control_flow| {
        if let Event::RedrawRequested(_) = event {
            let triangle_vertices: [Vector3D; 3] = [Vector3D::new(10., 10., 0.),
                Vector3D::new(100., 30., 0.), Vector3D::new(190., 160., 0.)];
            world.draw_triangle(pixels.get_frame(), &triangle_vertices);
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
            let y = (i / WIDTH as usize) as i16;
            
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
            // for x_index in (bbox_min.x as usize)..(bbox_max.x as usize) {
            //     for y_index in (bbox_min.y as usize)..(bbox_max.y as usize) {
            //         //     let x = (i % WIDTH as usize) as i16;
            //         //     let y = (i / WIDTH as usize) as i16;
            //         let temp_point = Vector3D::new(x_index as f32,
            //             y_index as f32, 0.0);
            //         let bc_screen = Vector3D::barycentric(pts, &temp_point);
            // 
            //         let outside_the_triangle: bool = bc_screen.x < 0.0
            //             || bc_screen.y < 0.0 || bc_screen.z < 0.0;
            //         
            //         let rgba = if outside_the_triangle {
            //             [0x5e, 0x48, 0xe8, 0xff]
            //         } else {
            //             [0x48, 0xb2, 0xe8, 0xff]
            //         };
            //         
            //         pixel.copy_from_slice(&rgba);
            //     }
            // }
        }
        // for (i, pixel) in frame.chunks_exact_mut(4).enumerate() {
        //     let x = (i % WIDTH as usize) as i16;
        //     let y = (i / WIDTH as usize) as i16;
        //     
        //     let inside_the_box = x >= self.box_x
        //         && x < self.box_x + BOX_SIZE
        //         && y >= self.box_y
        //         && y < self.box_y + BOX_SIZE;
        // 
        //     let rgba = if inside_the_box {
        //         [0x5e, 0x48, 0xe8, 0xff]
        //     } else {
        //         [0x48, 0xb2, 0xe8, 0xff]
        //     };
        // 
        //     pixel.copy_from_slice(&rgba);
        // }
    }
}