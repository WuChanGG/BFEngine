use crate::math::Vector3D;

mod math;

fn main() {
    let i = Vector3D::new(1.0, 0.0, 0.0);
    let j = Vector3D::new(0.0, 1.0, 0.0);
    let k = Vector3D::new(0.0, 0.0, 1.0);
    
    // let temp_cross = k.cross(&i);
    let temp_scalar_triple = i.cross(&j).dot(&k);
    println!("{:?}", temp_scalar_triple);
}
