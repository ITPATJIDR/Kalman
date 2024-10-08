import tkinter as tk
import serial
import math

# Initialize serial
port_name = '/dev/cu.usbmodem153480801'  # Update this with the appropriate port name
teensy = serial.Serial(port_name, 115200)

# Constants
triangle_radius = 30
scaling_factor = 10
line_buffer = ""

# Initialize Tkinter
root = tk.Tk()
root.attributes("-fullscreen", True)
canvas = tk.Canvas(root, width=1920, height=1080, bg='black')
canvas.pack()

def draw_led(x, y):
    x_scaled = x * scaling_factor
    y_scaled = y * scaling_factor
    canvas.create_rectangle(
        x_scaled - 10, y_scaled - 10, x_scaled + 10, y_scaled + 10,
        fill='white', outline='gray'
    )

def draw_device(x, y, theta):
    x1 = x * scaling_factor + triangle_radius * math.cos(math.pi / 3 + theta)
    y1 = y * scaling_factor + triangle_radius * math.sin(math.pi / 3 + theta)
    x2 = x * scaling_factor + triangle_radius * math.cos(math.pi + theta)
    y2 = y * scaling_factor + triangle_radius * math.sin(math.pi + theta)
    x3 = x * scaling_factor + triangle_radius * math.cos(-math.pi / 3 + theta)
    y3 = y * scaling_factor + triangle_radius * math.sin(-math.pi / 3 + theta)

    # Draw the triangle
    canvas.create_polygon(
        x1, y1, x2, y2, x3, y3,
        fill='red', outline='gray'
    )
    
    # Draw the line
    canvas.create_line(x1, y1, x3, y3, fill='green')

def update_canvas():
    global line_buffer

    if teensy.in_waiting > 0:
        char = teensy.read().decode()
        if char != '\n':
            line_buffer += char
        else:
            try:
                list_values = line_buffer.split()
                x = float(list_values[0])
                y = float(list_values[1])
                theta = float(list_values[2])

                # Clear canvas and redraw
                canvas.delete('all')
                draw_led(0, 0)
                draw_led(41, 0)
                draw_led(0, 41)
                draw_led(41, 41)
                draw_device(x * 100, y * 100, theta)

                line_buffer = ""
                teensy.reset_input_buffer()

            except (ValueError, IndexError):
                print("Invalid data received")

    # Schedule the next update
    root.after(8, update_canvas)  # Runs the function again after 8 milliseconds (~120 FPS)

# Exit on 'Escape' key press
def on_key(event):
    if event.keysym == 'Escape':
        root.destroy()

root.bind('<Key>', on_key)

# Start the update loop
update_canvas()
root.mainloop()
