use std::default;
use std::ops;
use bytemuck::*;

// A 3D math library made by me based on the "Foundations of game engine
// development: Mathemathics" book by Eric Lengyel

// The "image document" has all the images references here.
// The link to the image document is https://docs.google.com/document/d/11mnLJtuahx7EVwl0FdKtZ2hqah-MSPtywuxfZaUqZUI/edit

#[derive(Clone, Copy, Debug)]
#[repr(C)]
pub struct Vector3D
{
    pub x: f32,
    pub y: f32,
    pub z: f32,
}

unsafe impl Zeroable for Vector3D {}

unsafe impl Zeroable for &Vector3D {}

unsafe impl bytemuck::Pod for Vector3D {}

impl Vector3D {
    pub fn new(in_x: f32, in_y: f32, in_z: f32) -> Self {
        return Self {
            x: in_x,
            y: in_y,
            z: in_z,
        };
    }

    pub fn dot(self, rhs: &Vector3D) -> f32 {
        return self.x * rhs.x + self.y * rhs.y + self.z * rhs.z;
    }

    // When two vectors are parallel the cross product is the zero vector
    // When two vectors are not parallel the cross product is a vector that is
    // perpendicular to both "a" and "b"
    // ||axb||=||a||||b||sin(A) where A is the planar angle between
    // the directions of "a" and "b"
    // the area of the parallelogram having sides "a" and "b" is defined
    // by Area = ||axb||
    pub fn cross(self, rhs: &Vector3D) -> Vector3D {
        return Vector3D::new(self.y * rhs.z - self.z * rhs.y,
            self.z * rhs.x - self.x * rhs.z,
            self.x * rhs.y - self.y * rhs.x);
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

    fn project(self, other: &Vector3D) -> Vector3D {
        return *other * (self.dot(other) / other.dot(other));
    }

    fn reject(self, other: &Vector3D) -> Vector3D {
        return self - *other * (self.dot(other) / other.dot(other));
    }

    // u_i = v_i - v_i.project(u_k).summation()
    // == v_i - ( (v_i.dot(u_k) / u_k.squared) * u_k ).summation()
    // for example a set of three vectors {v1, v2, v3} is orthogonalized
    // by using the calculations
    // u_1 = v_1
    // u_2 = v_2 - v_2.project(u_1)
    // u_3 = v_3 - v_3.project(u_1) - v_3.project(u_2)
    // Vector u_i must be renormalized to unit length after orthogonalization
    // by dividing each one by its magnitude, this is called orthonormalization
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

impl ops::Mul<f32> for &Vector3D {
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

// The inverse matrix is only defined for square matrices, and the inverse
// of an nxn matrix is a matrix denoted by M^-1 having the property
// M * M^-1 = I_n (Identity matrix)

// Topic: Elementary Matrix, Elementary Matrices
// Elementary row operations, they affect one or two entire rows of a matrix
// These are the three operations
// - Multiply one row of M by a nonzero scalar value
// - Exchange two rows of M
// - Add a scalar multiple of one row of M to another row of M
// Each elementary row operation can be applied to a nxn matrix

// Application of rule 1:
// To multiply row r by a scalar value t, the ele
// mentary matrix E has the following form, where the (r,r) entry of the identity
// matrix has been replaced by t.

// Application of rule 2:
// To exchange row r and row s, the elementary matrix E has the following form,
// where the same rows have been exchanged in the identity matrix.

// Application of rule 3:
// To add row s multiplied by the scalar value t to row r, the elementary matrix E has
// the following forro, where the ( r, s) entry of the identity matrix has been replaced
// by t.

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

    // Topic: Determinant
    // The determinant of a nxn matrix M is a scalar value that can be thought as
    // the magnitude for M it is written "det(M)" or "|M|"
    // If we consider the n-columns or the n-rows of the matrix as a set of
    // vectors then the determinant is equal to the hypervolume of the n-dimensional
    // parallelotope formed by those vectors, and it can be positive or negative
    // A matrix has a determinant only if the determinant is non-zero

    // A term would be the collection {0,1,2} and each part of the collection is
    // a factor. A set of terms based on all possible permutations of {0,1,2}
    // indicate the columns (j index), and the row index is always {0,1,2}
    // ../images/MatrixDeterminantFormula.png
    fn determinant(self) -> f32 {
        return self.get_entry(0, 0) * (self.get_entry(1, 1) * self.get_entry(2, 2)
            - self.get_entry(1, 2) * self.get_entry(2, 1))
            + self.get_entry(0, 1) * (self.get_entry(1, 2) * self.get_entry(2, 0)
            - self.get_entry(1, 0) * self.get_entry(2, 2))
            + self.get_entry(0, 2) * (self.get_entry(1, 0) * self.get_entry(2, 1)
            - self.get_entry(1, 1) * self.get_entry(2, 0));
        // ax * (by * 1 - cy * 1) ==  
    }

    // Subtopic: Expansion by minors
    // The notation M_i_j (with a line above the indices) represent a submatrix
    // of M that excludes row "i" and column "j", the overbar (the line above)
    // is interpreted as "not"
    // ../images/determinant_formula_by_minors.png for more details

    fn inverse(&self) -> Matrix3D {
        let a: &Vector3D = self.get_vector_ref(0);
        let b: &Vector3D = self.get_vector_ref(1);
        let c: &Vector3D = self.get_vector_ref(2);

        let r0: Vector3D = b.cross(c);
        let r1: Vector3D = c.cross(a);
        let r2: Vector3D = a.cross(b);

        let inv_det: f32 = 1.0 / r2.dot(c);

        Matrix3D::new_from_vectors(
            &Vector3D::new(r0.x * inv_det, r0.y * inv_det, r0.z * inv_det),
            &Vector3D::new(r1.x * inv_det, r1.y * inv_det, r1.z * inv_det),
            &Vector3D::new(r2.x * inv_det, r2.y * inv_det, r2.z * inv_det),
        )
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

#[derive(Clone, Copy, Debug)]
#[repr(C)]
struct Matrix4D {
    entries: [[f32; 4]; 4],
}

impl Matrix4D {
    fn inverse(&self) -> Matrix4D {
        let a: &Vector3D = try_cast_ref::<[f32; 4], Vector3D>(
            &self.entries[0]).unwrap();
        let b: &Vector3D = try_cast_ref::<[f32; 4], Vector3D>(
            &self.entries[1]).unwrap();
        let c: &Vector3D = try_cast_ref::<[f32; 4], Vector3D>(
            &self.entries[2]).unwrap();
        let d: &Vector3D = try_cast_ref::<[f32; 4], Vector3D>(
            &self.entries[3]).unwrap();

        let x: f32 = self.entries[3][0];
        let y: f32 = self.entries[3][1];
        let z: f32 = self.entries[3][2];
        let w: f32 = self.entries[3][3];

        let mut s: Vector3D = a.cross(b);
        let mut t: Vector3D = c.cross(d);
        let mut u: Vector3D = a * y - b * x;
        let mut v: Vector3D = c * w - d * z;

        let inv_det: f32 = 1.0 / s.dot(&v) + t.dot(&u);
        s *= inv_det;
        t *= inv_det;
        u *= inv_det;
        v *= inv_det;

        let r0: Vector3D = b.cross(&v) + t * y;
        let r1: Vector3D = v.cross(&a) - t * x;
        let r2: Vector3D = d.cross(&u) + s * w;
        let r3: Vector3D = u.cross(&c) - s * z;

        return Matrix4D::new(
            r0.x, r0.y, r0.z, -b.dot(&t),
            r1.x, r1.y, r1.z, a.dot(&t),
            r2.x, r2.y, r2.z, -d.dot(&s),
            r3.x, r3.y, r3.z, c.dot(&s),
        );
    }

    fn new(n00: f32, n01: f32, n02: f32, n03: f32,
        n10: f32, n11: f32, n12: f32, n13: f32,
        n20: f32, n21: f32, n22: f32, n23: f32,
        n30: f32, n31: f32, n32: f32, n33: f32,
    ) -> Self {
        let mut zero_matrix: Matrix4D = Matrix4D::default();
        zero_matrix.entries[0][0] = n00;
        zero_matrix.entries[0][1] = n10;
        zero_matrix.entries[0][2] = n20;
        zero_matrix.entries[0][3] = n30;
        zero_matrix.entries[1][0] = n01;
        zero_matrix.entries[1][1] = n11;
        zero_matrix.entries[1][2] = n21;
        zero_matrix.entries[1][3] = n31;
        zero_matrix.entries[2][0] = n02;
        zero_matrix.entries[2][1] = n12;
        zero_matrix.entries[2][2] = n22;
        zero_matrix.entries[2][3] = n32;
        zero_matrix.entries[3][0] = n03;
        zero_matrix.entries[3][1] = n13;
        zero_matrix.entries[3][2] = n23;
        zero_matrix.entries[3][3] = n33;
        return zero_matrix;
    }
}

impl Default for Matrix4D {
    fn default() -> Self {
        let temp_entries: [[f32; 4]; 4] = [[0.0; 4]; 4];
        Matrix4D { entries: temp_entries }
    }
}

impl ops::BitXor<Vector3D> for Vector3D {
    type Output = Vector3D;
    fn bitxor(self, rhs: Vector3D) -> Self::Output {
        Vector3D::new(self.x.powf(rhs.x), self.y.powf(rhs.y),
            self.z.powf(rhs.z))
    }
}

static WIDTH: i32 = 200;
static HEIGHT: i32 = 200;

impl Vector3D {
    // compute 2D determinant
    fn orient_2d(a: &Vector3D, b: &Vector3D, c: &Vector3D) -> f32 {
        return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
    }

    fn draw_triagle(v0: &Vector3D, v1: &Vector3D, v2: &Vector3D) {
        // compute the triangle bounding box
        let mut min_x = v0.x.min(v1.x).min(v2.x);
        let mut min_y = v0.y.min(v1.y).min(v2.y);
        let mut max_x = v0.x.max(v1.x).max(v2.x);
        let mut max_y = v0.y.max(v1.y).max(v2.y);

        // clip against screen bounds
        min_x = min_x.max(0.0);
        min_y = min_y.max(0.0);
        max_x = max_x.min(WIDTH as f32 - 1);
        max_y = max_y.min(HEIGHT as f32 - 1);

        let min_y_int: i32 = min_y as i32;
        let max_y_int: i32 = max_y as i32;
        let min_x_int: i32 = min_x as i32;
        let max_x_int: i32 = max_x as i32;

        // Rasterize
        for point_y in min_y_int..max_y_int {
            for point_x in min_x_int..max_x_int {
                let p: Vector3D = Vector3D::new(point_x as f32,
                    point_y as f32, 0.0);
                
                // determine barycentric coordinate
                let w0: f32 = orient2d(v1, v2, p);
                let w1: f32 = orient2d(v2, v0, p);
                let w2: f32 = orient2d(v0, v1, p);
                
                if (w0 >= 0. && w1 >= 0. && w2 >= 0.)
                {
                    // Render pixel
                }
            }
        }
    }
    
    // https://users.csc.calpoly.edu/~zwood/teaching/csc471/2017F/barycentric.pdf
    fn barycentric(pts: &[Vector3D; 3], p: &Vector3D) {
        let area_of_triangle: f32 = (pts[1][0] - pts[0][0]) * (pts[2][1] - pts[0][1])
            - (pts[2][0] - pts[0][0]) * (pts[1][1] - pts[0][1]);
        
        let u_nominator: f32 = (pts[0][0] - pts[2][0]) * (p[1] - pts[2][1])
            - (p[0] - pts[2][0]) * (pts[0][1] - pts[2][1]);
        
        // the "weight" u, β = u in the formula below
        // (1 − β − γ) a + βb + γc
        // same applies for the weight "w", γ = w
        let u_nominator: f32 = (pts[0][0] - pts[2][0]) * (p[1] - pts[2][1])
            - (p[0] - pts[2][0]) * (pts[0][1] - pts[2][1]);
        
        let v_nominator: f32 = (pts[0][0] - pts[2][0]) * (p[1] - pts[2][1])
            - (p[0] - pts[2][0]) * (pts[0][1] - pts[2][1]);
        
        let u: f32 = u_nominator / area_of_triangle;
    }
}

// Vec3f barycentric(Vec2i *pts, Vec2i P) {
// Vec3f u = Vec3f(pts[2][0]-pts[0][0], pts[1][0]-pts[0][0], pts[0][0]-P[0])
//  ^Vec3f(pts[2][1]-pts[0][1], pts[1][1]-pts[0][1], pts[0][1]-P[1]);
// /* `pts` and `P` has integer value as coordinates
//    so `abs(u[2])` < 1 means `u[2]` is 0, that means
//    triangle is degenerate, in this case return something with negative coordinates */
// if (std::abs(u.z)<1) return Vec3f(-1,1,1);
// return Vec3f(1.f-(u.x+u.y)/u.z, u.y/u.z, u.x/u.z);
// } 