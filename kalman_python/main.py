import pygame
import serial
import math

# Initialize serial
port_name = '/dev/cu.usbmodem153480801'  # Update this with the appropriate port name
teensy = serial.Serial(port_name, 115200)

# Pygame setup
pygame.init()
screen = pygame.display.set_mode((1920, 1080), pygame.RESIZABLE)
clock = pygame.time.Clock()
font = pygame.font.Font(None, 36)

# Constants
triangle_radius = 30
scaling_factor = 10
line_buffer = ""

def draw_led(x, y):
    rect = pygame.Rect(x * scaling_factor - 10, y * scaling_factor - 10, 20, 20)
    pygame.draw.rect(screen, (255, 255, 255), rect)
    pygame.draw.rect(screen, (128, 128, 128), rect, 1)

def draw_device(x, y, theta):

    y = -y
    y += 30
    x1 = x * scaling_factor + triangle_radius * math.cos(math.pi + (math.pi / 3) + theta)
    y1 = y * scaling_factor + triangle_radius * math.sin(math.pi + (math.pi / 3) + theta)
    x2 = x * scaling_factor + triangle_radius * math.cos(math.pi + theta)
    y2 = y * scaling_factor + triangle_radius * math.sin(math.pi + theta)
    x3 = x * scaling_factor + triangle_radius * math.cos(math.pi + (-math.pi / 3) + theta)
    y3 = y * scaling_factor + triangle_radius * math.sin(math.pi + (-math.pi / 3) + theta)
    
    pygame.draw.polygon(screen, (255, 0, 0), [(x1, y1), (x2, y2), (x3, y3)])
    pygame.draw.line(screen, (0, 255, 0), (x1, y1), (x3, y3), 1)

# Main loop
running = True
while running:
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False
        elif event.type == pygame.KEYDOWN and event.key == pygame.K_ESCAPE:
            running = False

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

                screen.fill((0, 0, 0))
                screen_center = (screen.get_width() // 4, screen.get_height() // 2)
                pygame.draw.line(screen, (255, 255, 255), (0, screen_center[1]), (screen.get_width(), screen_center[1]), 1)

                draw_led(0, 0)
                draw_led(41, 0)
                draw_led(0, 41)
                draw_led(41, 41)
                draw_device(x * 100, y * 100, theta)

                pygame.display.flip()
                line_buffer = ""
                teensy.reset_input_buffer()

            except (ValueError, IndexError):
                print("Invalid data received")

    clock.tick(120)

pygame.quit()
