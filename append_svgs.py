import argparse
import xml.etree.ElementTree as ET

def get_svg_dimensions(root):
    viewBox = root.attrib.get('viewBox')
    if viewBox:
        return list(map(float, viewBox.split()))
    else:
        width = float(root.attrib.get('width', '0').rstrip('px'))
        height = float(root.attrib.get('height', '0').rstrip('px'))
        return [0, 0, width, height]

def append_svg_legend(main_file, legend_file, gap=10, output_file=None):
    # Parse SVG files
    ET.register_namespace("", "http://www.w3.org/2000/svg")
    main_tree = ET.parse(main_file)
    main_root = main_tree.getroot()
    legend_tree = ET.parse(legend_file)
    legend_root = legend_tree.getroot()

    # Get dimensions
    main_x, main_y, main_width, main_height = get_svg_dimensions(main_root)
    legend_x, legend_y, legend_width, legend_height = get_svg_dimensions(legend_root)

    # Calculate new viewBox
    new_width = max(main_width, legend_width)
    new_height = main_height + gap + legend_height
    new_viewBox = f"{main_x} {main_y} {new_width} {new_height}"

    # Update main SVG dimensions
    main_root.set('viewBox', new_viewBox)
    main_root.set('width', str(new_width))
    main_root.set('height', str(new_height))

    # Create a group for the legend
    legend_group = ET.SubElement(main_root, 'g')
    legend_x_pos = main_x + main_width - legend_width  # Right-align
    legend_y_pos = main_y + main_height + gap
    legend_group.set('transform', f'translate({legend_x_pos}, {legend_y_pos})')

    # Copy legend content to the new group
    for child in legend_root:
        legend_group.append(child)

    # Save the combined SVG
    if output_file is None:
        output_file = main_file
    tree = ET.ElementTree(main_root)
    tree.write(output_file, encoding='unicode', xml_declaration=True)
    print(f"Combined SVG saved as {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Combine a main SVG with a legend SVG.")
    parser.add_argument("main_file", help="Path to the main SVG file")
    parser.add_argument("legend_file", help="Path to the legend SVG file")
    parser.add_argument("--gap", type=int, default=10, help="Gap between main image and legend (default: 10)")
    parser.add_argument("--output", default=None, help="Output file name (default: overwrite main file)")

    args = parser.parse_args()
    append_svg_legend(args.main_file, args.legend_file, args.gap, args.output)

if __name__ == "__main__":
    main()
