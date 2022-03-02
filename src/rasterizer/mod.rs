use crate::math::{Vector3D};

static width: i32 = 200;
static height: i32 = 200;

impl std::vec {
    fn barycentric(pts: [Vector3D; 3], point: Vector3D) -> Vector3D {
        let u: Vector3D = Vector3D::new(
            pts[2][0] - pts[0][0], pts[1][0] - pts[0][0], pts[0][0] - point[0]
        )^Vector3D::new(
            pts[2][1] - pts[0][1], pts[1][1] - pts[0][1], pts[0][1] - point[1]);
        // Vec3f u = Vec3f(pts[2][0]-pts[0][0], pts[1][0]-pts[0][0],
        // pts[0][0]-P[0])^Vec3f(pts[2][1]-pts[0][1], pts[1][1]-pts[0][1], pts[0][1]-P[1]);
        return u;
    }
}