//#[path="common/common_run.rs"]
//mod wgpu_run01;
mod common;

use winit::{
    event_loop::{EventLoop},
};
use std::borrow::Cow;

fn main() {
    let event_loop = EventLoop::new();
    let window = winit::window::Window::new(&event_loop).unwrap(); 
    window.set_title("Piccolo Tirangle");
    env_logger::init();

    let inputs = common::Inputs{
        source: wgpu::ShaderSource::Wgsl(Cow::Borrowed(include_str!("../shader/first_triangle.wgsl"))),
        topology: wgpu::PrimitiveTopology::TriangleList,
        strip_index_format: None
    };

    pollster::block_on( common::run(event_loop, window, inputs, 3));    
}
