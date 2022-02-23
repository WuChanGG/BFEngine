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

impl ops::IndexMut<usize> for Vector3D {
    fn index_mut<'a>(&'a mut self, i: usize)
        -> &'a mut f32 {
        match i {
            0 => return &mut self.x,
            1 => return &mut self.y,
            2 => return &mut self.z,
            _ => panic!("Vector::Index index out of bounds")
        }
    }
}

impl ops::MulAssign<f32> for Vector3D {
    fn mul_assign(&mut self, rhs: f32) {
        self.x *= rhs;
        self.y *= rhs;
        self.z *= rhs;
    }
}

impl ops::DivAssign<f32> for Vector3D {
    fn div_assign(&mut self, mut rhs: f32) {
        rhs = 1.0 / rhs;
        self.x *= rhs;
        self.y *= rhs;
        self.z *= rhs;
    }
}

impl ops::Mul<f32> for Vector3D {
    type Output = Vector3D;
    fn mul(self, rhs: f32) -> Self::Output {
        return Vector3D {
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z *rhs,
        };
    }
}

impl ops::Div<f32> for &Vector3D {
    type Output = Vector3D;
    fn div(self, mut rhs: f32) -> Self::Output {
        rhs = 1.0 / rhs;
        return Vector3D {
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z *rhs,
        };
    }
}

impl ops::Div<f32> for Vector3D {
    type Output = Vector3D;
    fn div(self, mut rhs: f32) -> Self::Output {
        rhs = 1.0 / rhs;
        return Vector3D {
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z *rhs,
        };
    }
}

impl ops::Neg for Vector3D {
    type Output = Vector3D;
    fn neg(self) -> Self::Output {
        return Vector3D {
            x: -self.x,
            y: -self.y,
            z: -self.z
        };
    }
}

impl Vector3D {
    fn magnitude(&self) -> f32 {
        return f32::sqrt(
            self.x * self.x + self.y * self.y + self.z * self.z)
    }

    fn normalize(&self) -> Vector3D {
        return self / self.magnitude();
    }
}