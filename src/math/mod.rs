use std::default;
use std::ops;

// A 3D math library made by me based on the "Foundations of game engine
// development: Mathemathics" book by Eric Lengyel

// The "image document" has all the images references here.
// The link to the image document is https://docs.google.com/document/d/11mnLJtuahx7EVwl0FdKtZ2hqah-MSPtywuxfZaUqZUI/edit

struct Vector3D
{
    x : f32,
    y : f32,
    z : f32
}

impl Vector3D {
    fn new (in_x: f32, in_y: f32, in_z: f32) -> Self {
        return Self {
            x: in_x,
            y: in_y,
            z: in_z
        };
    }
}

impl Default for Vector3D {
    fn default () -> Self {
        return Self {
            x: 0.0,
            y: 0.0,
            z: 0.0
        }
    }
}

impl ops::Index<usize> for Vector3D {
    type Output = f32;
    fn index<'a>(&'a self, i: usize)
        -> &'a f32 {
        match i {
            0 => return &self.x,
            1 => return &self.y,
            2 => return &self.z,
            _ => panic!("Vector::Index index out of bounds")
        }
    }
}