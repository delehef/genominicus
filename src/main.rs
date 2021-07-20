use crate::nhx::*;

mod nhx;
fn main() {
    let r = Tree::from_filename("../done/prims_065/final_trees/76-profilenj.nhx").unwrap();
    r.print();
}
