use std::default;
use std::ops;
use bytemuck::*;

// A 3D math library made by me based on the "Foundations of game engine
// development: Mathemathics" book by Eric Lengyel

// The "image document" has all the images references here.
// The link to the image document is https://docs.google.com/document/d/11mnLJtuahx7EVwl0FdKtZ2hqah-MSPtywuxfZaUqZUI/edit

#[derive(Clone, Copy, Debug)]
#[repr(C)]
struct Vector3D
{
    x: f32,
    y: f32,
    z: f32,
}

unsafe impl Zeroable for Vector3D {}

unsafe impl Zeroable for &Vector3D {}

unsafe impl bytemuck::Pod for Vector3D {}

impl Vector3D {
    fn new(in_x: f32, in_y: f32, in_z: f32) -> Self {
        return Self {
            x: in_x,
            y: in_y,
            z: in_z,
        };
    }

    fn dot(self, rhs: Vector3D) -> f32 {
        return self.x * rhs.x + self.y * rhs.y + self.z * rhs.z;
    }
    
    // When two vectors are parallel the cross product is the zero vector
    // When two vectors are not parallel the cross product is a vector that is
    // perpendicular to both "a" and "b"
    // ||axb||=||a||||b||sin(A) where A is the planar angle between
    // the directions of "a" and "b"
    // the area of the parallelogram having sides "a" and "b" is defined
    // by Area = ||axb||
    fn cross(self, rhs: Vector3D) -> Vector3D {
        return Vector3D::new(self.y * rhs.z - self.z * rhs.y,
        self.z * rhs.x - self.x * rhs.z,
        self.x * rhs.y - self.z * rhs.y)
    }
    
    // Topic: Scalar Triple Product.
    // Scalar triple product. The notation is
    // [a,b,c] == a.cross(b).dot(c) == b.cross(c).dot(a) == c.cross(b).dot(a)
    // the order of the vectors do not matter, they can wrap, and if they are
    // wrapped in the opposite direction then the scalar triple product is
    // negated, and this accounts for all permutations:
    // [c,b,a] == c.cross(b).dot(a) == b.cross(a).dot(c) == a.cross(b).dot(c)
    // == -[a,b,c]
    
    // The scalar triple product a.cross(b).dot(c) yields the volume of the
    // parallelepied (a 3D parallelogram) spanned by vectors a,b,c
    
    // area of the base is
    // Area = ||axb||
    
    // Height = ||c||sin(AngleX) where the angle goes from the
    // plane formed by "a" and "b" until intercepting the vector "c"
    // if this angle is measured by the angle between "c" and the line
    // a.cross(b) then Height = ||c||cos(AngleY), AngleX and AngleY are
    // complementary angles
    
    // Topic: Vector Projection.
    // Given a particular vector, find two or more other vectors with
    // specific alignments that add up to our original vector "decomposing
    // a vector into its separate components"
    // Each coordinate is equal to the magnitude of vector "v" multiplied
    // by the cosine of the angle that "v" makes with the corresponding axis
    // this means --> v = v.dot(i)*i + v.dot(j)*j + v.dot(k)*k
    // Image can be seen on page 41
    
    // In general the dot product can be used to project any vector "a"
    // onto any other non-zero vector "b" with the formula
    // a_projected_to_b = ( (a.dot(b)) / b.squared() ) * b
    // this notation indicates the component of the vector "a" that is 
    // parallel to the vector "b" and the equation gives us projection of
    // "a" onto "b"
    // Also if a.dot(b) (in the formula) is negative the a.projected_to(b)
    // is still parallel but in the opposite direction
    
    // the projection of a onto b can be expressed as the matrix product
    // (1 / vector_b.squared) * matrix_b * matrix_b_transposed * vector_a
    // this equation is an example of an outer product
    // the outer product between two vectors "u" and "v"
    
    // if we substract the projection a.projected_to(b) from the original
    // vector "a" then we get the part that is perpendicular to vector "b"
    // because we removed everything that is parallel to "b"
    // the perpendicular part of the decomposition is called rejection of "a"
    // from "b" and is written
    // a.rejection_from(b) = a - a.projected_to(b)
    // = a - ( (a.dot(b) / b.squared) ) * b
    
    // a.projected_to(b) and a.rejected_from(b) form the sides of a right
    // triangle where "a" is the hypotenuse so
    // a.projected_to(b).squared + a.rejected_from(b).squared = a.squared
    // and basic geometry tells us that
    // a.projected_to(b).absolute = a.absolute().cos(angle)
    // a.rejected_from(b).absolute() = a.absolute().sin(angle)
    
    // Applications of vector projection:
    // Orthogonalization: in which each member in a set of vectors is modified
    // so that it is perpendicular, or orthogonal to all other vectors.
    // E.G. if we have two vector "a" and "b": 
    // replacing "a" with a.rejected_from(b) or b.rejected_from(a) for "b"
    // we are substracting the projection of one vector onto the other so that
    // the parallel component is removed leaving only the perpendicular
    // component
    
    fn project(self, other: Vector3D) -> Vector3D {
        return other * (self.dot(other) / other.dot(other));
    }
    
    fn reject(self, other: Vector3D) -> Vector3D {
        return self - other * (self.dot(other) / other.dot(other));
    }
    
    // u_i = v_i - v_i.project(u_k).summation()
}

impl Default for Vector3D {
    fn default() -> Self {
        return Self {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        };
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
            z: self.z * rhs,
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
            z: self.z * rhs,
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
            z: self.z * rhs,
        };
    }
}

impl ops::Neg for Vector3D {
    type Output = Vector3D;
    fn neg(self) -> Self::Output {
        return Vector3D {
            x: -self.x,
            y: -self.y,
            z: -self.z,
        };
    }
}

impl Vector3D {
    fn magnitude(&self) -> f32 {
        return f32::sqrt(
            self.x * self.x + self.y * self.y + self.z * self.z);
    }

    fn normalize(&self) -> Vector3D {
        return self / self.magnitude();
    }
}

impl ops::AddAssign<Vector3D> for Vector3D {
    fn add_assign(&mut self, rhs: Vector3D) {
        self.x += rhs.x;
        self.y += rhs.y;
        self.z += rhs.z;
    }
}

impl ops::SubAssign<Vector3D> for Vector3D {
    fn sub_assign(&mut self, rhs: Vector3D) {
        self.x -= rhs.x;
        self.y -= rhs.y;
        self.z -= rhs.z;
    }
}

impl ops::Add<Vector3D> for Vector3D {
    type Output = Vector3D;
    fn add(self, rhs: Vector3D) -> Vector3D {
        Vector3D::new(self.x + rhs.x, self.y + rhs.y,
            self.z + rhs.z)
    }
}

impl ops::Sub<Vector3D> for Vector3D {
    type Output = Vector3D;
    fn sub(self, rhs: Vector3D) -> Vector3D {
        Vector3D::new(self.x - rhs.x, self.y - rhs.y,
            self.z - rhs.z)
    }
}

// Matrix Section
// - The numbers that make up a matrix M are called entries

// - When every entry is 0 except the diagonal entries, the matrix is called
// "diagonal"

// - A "transpose" of a matrix "M" with size "n x m" is denoted by "M^T"
// and has size "m x n"

// - If a matrix is equal to the transpose, such as M_i_j == M_j_i
// the matrix is called symmetric

// - If the entries of M^t are equal to the negations of the same entries
// in the matrix M (M_i_j^T == -Mij) then the matrix M is called antisymmetric.
// for this to be the case all diagonal entries must be zero

struct Matrix3D {
    entries: [[f32; 3]; 3],
}

impl Matrix3D {
    fn new(n00: f32, n01: f32, n02: f32,
        n10: f32, n11: f32, n12: f32,
        n20: f32, n21: f32, n22: f32,
    ) -> Self {
        let mut zero_matrix: Matrix3D = Matrix3D::default();
        zero_matrix.entries[0][0] = n00;
        zero_matrix.entries[0][1] = n10;
        zero_matrix.entries[0][2] = n20;
        zero_matrix.entries[1][0] = n01;
        zero_matrix.entries[1][1] = n11;
        zero_matrix.entries[1][2] = n21;
        zero_matrix.entries[2][0] = n02;
        zero_matrix.entries[2][1] = n12;
        zero_matrix.entries[2][2] = n22;
        return zero_matrix;
    }

    fn new_from_vectors(a: &Vector3D, b: &Vector3D, c: &Vector3D) -> Self {
        let mut temp_matrix: Matrix3D = Matrix3D::default();
        temp_matrix.entries[0][0] = a.x;
        temp_matrix.entries[0][1] = a.y;
        temp_matrix.entries[0][2] = a.z;
        temp_matrix.entries[1][0] = b.x;
        temp_matrix.entries[1][1] = b.y;
        temp_matrix.entries[1][2] = b.z;
        temp_matrix.entries[2][0] = c.x;
        temp_matrix.entries[2][1] = c.y;
        temp_matrix.entries[2][2] = c.z;
        return temp_matrix;
    }

    fn get_entry(&self, i: usize, j: usize) -> f32 {
        return self.entries[j][i];
    }

    fn get_mut_entry(&mut self, i: usize, j: usize) -> &mut f32 {
        return &mut self.entries[j][i];
    }

    fn get_ref_entry(&self, i: usize, j: usize) -> &f32 {
        return &self.entries[j][i];
    }

    fn get_vector_mut(&mut self, j: usize) -> &mut Vector3D {
        let mut temp: &mut Vector3D = try_cast_mut::<[f32; 3], Vector3D>(
            &mut self.entries[j]).unwrap();
        return temp;
    }

    fn get_vector_ref(&self, j: usize) -> &Vector3D {
        let temp: &Vector3D = try_cast_ref::<[f32; 3], Vector3D>(
            &self.entries[j]).unwrap();
        return temp;
    }
}

impl Default for Matrix3D {
    fn default() -> Self {
        let temp_entries: [[f32; 3]; 3] = [[0.0; 3]; 3];
        Matrix3D { entries: temp_entries }
    }
}

impl ops::Mul<Matrix3D> for Matrix3D {
    type Output = Matrix3D;
    fn mul(self, rhs: Matrix3D) -> Matrix3D {
        return Matrix3D::new(
            // n00
            self.get_entry(0, 0) * rhs.get_entry(0, 0)
                + self.get_entry(0, 1) * rhs.get_entry(1, 0)
                + self.get_entry(0, 2) * rhs.get_entry(2, 0),
            // n01
            self.get_entry(0, 0) * rhs.get_entry(0, 1)
                + self.get_entry(0, 1) * rhs.get_entry(1, 1)
                + self.get_entry(0, 2) * rhs.get_entry(2, 1),
            // n02
            self.get_entry(0, 0) * rhs.get_entry(0, 2)
                + self.get_entry(0, 1) * rhs.get_entry(1, 2)
                + self.get_entry(0, 2) * rhs.get_entry(2, 2),
            // n10
            self.get_entry(1, 0) * rhs.get_entry(0, 0)
                + self.get_entry(1, 1) * rhs.get_entry(1, 0)
                + self.get_entry(1, 2) * rhs.get_entry(2, 0),
            // n11
            self.get_entry(1, 0) * rhs.get_entry(0, 1)
                + self.get_entry(1, 1) * rhs.get_entry(1, 1)
                + self.get_entry(1, 2) * rhs.get_entry(2, 1),
            // n12
            self.get_entry(1, 0) * rhs.get_entry(0, 2)
                + self.get_entry(1, 1) * rhs.get_entry(1, 2)
                + self.get_entry(1, 2) * rhs.get_entry(2, 2),
            // n20
            self.get_entry(2, 0) * rhs.get_entry(0, 0)
                + self.get_entry(2, 1) * rhs.get_entry(1, 0)
                + self.get_entry(2, 2) * rhs.get_entry(2, 0),
            // n21
            self.get_entry(2, 0) * rhs.get_entry(0, 1)
                + self.get_entry(2, 1) * rhs.get_entry(1, 1)
                + self.get_entry(2, 2) * rhs.get_entry(2, 1),
            // n22
            self.get_entry(2, 0) * rhs.get_entry(0, 2)
                + self.get_entry(2, 1) * rhs.get_entry(1, 2)
                + self.get_entry(2, 2) * rhs.get_entry(2, 2),
        );
    }
}

impl ops::Mul<Vector3D> for Matrix3D {
    type Output = Vector3D;
    fn mul(self, rhs: Vector3D) -> Vector3D {
        return Vector3D::new(
            self.get_entry(0, 0) * rhs.x
                + self.get_entry(0, 1) * rhs.y
                + self.get_entry(0, 2) * rhs.z,
            self.get_entry(1, 0) * rhs.x
                + self.get_entry(1, 1) * rhs.y
                + self.get_entry(1, 2) * rhs.z,
            self.get_entry(1, 0) * rhs.x
                + self.get_entry(1, 1) * rhs.y
                + self.get_entry(1, 2) * rhs.z,
        );
    }
}
