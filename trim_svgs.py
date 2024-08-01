#!/usr/bin/env python
import os
import sys
import xml.etree.ElementTree as ET
from PIL import ImageFont, ImageDraw, Image

def get_float(element, attribute, default=0.0):
    return float(element.get(attribute, default))

def get_text_bbox(text, font_size=16, font_family='Verdana'):
    """
    Get the bounding box of a text element using Pillow's ImageDraw.textbbox()
    
    :param text: The text content
    :param font_size: Font size in pixels
    :param font_family: Font family name
    :return: A tuple (left, top, right, bottom) in pixels
    """
    try:
        font = ImageFont.truetype(f"{font_family}.ttf", font_size)
    except IOError:
        print(f"Warning: Font {font_family} not found. Using default font.")
        font = ImageFont.load_default()
    
    # Create a dummy image (we don't need an actual image)
    img = Image.new('RGB', (1, 1))
    draw = ImageDraw.Draw(img)
    
    # Get the bounding box
    bbox = draw.textbbox((0, 0), text, font=font)
    return bbox

def trim_svg_file(svg_file):
    # Parse the SVG file
    tree = ET.parse(svg_file)
    root = tree.getroot()

    # Create a new group element to hold the highlight rectangles
    highlight_group = ET.Element('g', {'id': 'text-highlights'})
    root.insert(0, highlight_group)  # Insert at the beginning to be behind other elements

    # Get the bounding box of the SVG content
    min_x = float('inf')
    min_y = float('inf')
    max_x = float('-inf')
    max_y = float('-inf')

    for element in root.iter():
        if element.tag.endswith(('rect', 'circle', 'ellipse', 'line', 'polyline', 'polygon', 'path', 'text')):
            if element.tag.endswith('rect'):
                x = get_float(element, 'x')
                y = get_float(element, 'y')
                width = get_float(element, 'width')
                height = get_float(element, 'height')
                min_x = min(min_x, x)
                min_y = min(min_y, y)
                max_x = max(max_x, x + width)
                max_y = max(max_y, y + height)
            elif element.tag.endswith('text'):
                x = get_float(element, 'x')
                y = get_float(element, 'y')
                font_size = get_float(element, 'font-size', 16)
                font_family = element.get('font-family', 'Verdana')
                text = element.text or ''
                anchor = element.get('text-anchor', 'start')
                
                left, top, right, bottom = get_text_bbox(text, font_size, font_family)
                width = right - left
                height = bottom - top

                width *= 1.5
                
                # Adjust x based on text-anchor
                if anchor == 'middle':
                    x -= width / 2
                elif anchor == 'end':
                    x -= width

                print(f"Text: '{text}', X: {x}, Y: {y}, Width: {width}, Height: {height}")  # Debug output

                # # Create a rectangle element
                # rect = ET.Element('rect', {
                #     'x': str(x),
                #     'y': str(y - height),  # Adjust y to align with text baseline
                #     'width': str(width),
                #     'height': str(height),
                #     'fill': 'red',
                #     'opacity': '0.3'  # Make it semi-transparent
                # })
                # highlight_group.append(rect)
                
                min_x = min(min_x, x)
                min_y = min(min_y, y - height)
                max_x = max(max_x, x + width)
                max_y = max(max_y, y)
            elif element.tag.endswith('circle'):
                cx = get_float(element, 'cx')
                cy = get_float(element, 'cy')
                r = get_float(element, 'r')
                min_x = min(min_x, cx - r)
                min_y = min(min_y, cy - r)
                max_x = max(max_x, cx + r)
                max_y = max(max_y, cy + r)
            elif element.tag.endswith('ellipse'):
                cx = get_float(element, 'cx')
                cy = get_float(element, 'cy')
                rx = get_float(element, 'rx')
                ry = get_float(element, 'ry')
                min_x = min(min_x, cx - rx)
                min_y = min(min_y, cy - ry)
                max_x = max(max_x, cx + rx)
                max_y = max(max_y, cy + ry)
            elif element.tag.endswith('line'):
                x1 = get_float(element, 'x1')
                y1 = get_float(element, 'y1')
                x2 = get_float(element, 'x2')
                y2 = get_float(element, 'y2')
                min_x = min(min_x, x1, x2)
                min_y = min(min_y, y1, y2)
                max_x = max(max_x, x1, x2)
                max_y = max(max_y, y1, y2)
            elif element.tag.endswith(('polyline', 'polygon')):
                points = element.get('points')
                if points:
                    coordinates = points.split()
                    for i in range(0, len(coordinates), 2):
                        x = float(coordinates[i])
                        y = float(coordinates[i + 1])
                        min_x = min(min_x, x)
                        min_y = min(min_y, y)
                        max_x = max(max_x, x)
                        max_y = max(max_y, y)
            elif element.tag.endswith('path'):
                # Parsing path data is more complex and may require a separate library
                # For simplicity, let's assume the path's bounding box is defined by the 'd' attribute
                d = element.get('d')
                if d:
                    path_data = d.split()
                    i = 0
                    while i < len(path_data):
                        if path_data[i] in ['M', 'L']:
                            x = float(path_data[i + 1])
                            y = float(path_data[i + 2])
                            min_x = min(min_x, x)
                            min_y = min(min_y, y)
                            max_x = max(max_x, x)
                            max_y = max(max_y, y)
                            i += 3
                        else:
                            i += 1

    print(f"Calculated bounding box: ({min_x}, {min_y}, {max_x}, {max_y})")  # Debug output

    # Calculate new dimensions with padding
    padding = 5
    width = (max_x - min_x) + (2 * padding)
    height = (max_y - min_y) + (2 * padding)
    viewbox = f"{min_x - padding} {min_y - padding} {width} {height}"
    root.set('viewBox', viewbox)

    # Update the width and height attributes
    root.set('width', str(width))
    root.set('height', str(height))

    # Save the updated SVG file in place
    tree.write(svg_file)
    print(f"Updated viewBox: {viewbox}")
    print(f"Updated SVG: {svg_file}")

def main():
    # Check if filenames are provided as arguments
    if len(sys.argv) < 2:
        print("Please provide one or more SVG filenames as arguments.")
        sys.exit(1)

    # Get the list of SVG filenames from the command-line arguments
    svg_files = sys.argv[1:]

    # Process each SVG file
    for svg_file in svg_files:
        if not os.path.isfile(svg_file):
            print(f"File not found: {svg_file}")
            continue
        trim_svg_file(svg_file)

if __name__ == "__main__":
    main()
