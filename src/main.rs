use winit::{
    event::{Event, WindowEvent},
    event_loop::{ControlFlow, EventLoop},
    window::Window,
};

fn main() {
    println!("Hello, Piccolo!");

    let event_loop = EventLoop::new();
    let window = Window::new(&event_loop).unwrap();
    window.set_title("Piccolo Window");

    env_logger::init();

    // main loop
    event_loop.run(move |event, _, control_flow| {
        *control_flow = ControlFlow::Wait;
        match event {
            Event::WindowEvent {
                event: WindowEvent::CloseRequested,
                ..
            } => *control_flow = ControlFlow::Exit,
            _ => {}
        }
    });
}
