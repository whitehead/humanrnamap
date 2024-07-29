#!/usr/bin/env python
import argparse
from svgutils.compose import Figure, SVG
from svgutils.transform import fromfile

def append_svg_legend(main_file, legend_file, vertical='bottom', horizontal='right', gap=10, output_file=None):
    # Load the main visualization and the legend
    main_svg = fromfile(main_file)
    legend_svg = fromfile(legend_file)
    if output_file is None:
        output_file = main_file

    # Get the dimensions
    main_width = float(main_svg.width)
    main_height = float(main_svg.height)
    legend_width = float(legend_svg.width)
    legend_height = float(legend_svg.height)

    # Calculate the new dimensions and positions
    if vertical == 'bottom':
        new_height = main_height + legend_height
        main_y = 0
        legend_y = main_height + gap
    else:  # top
        new_height = main_height + legend_height
        main_y = legend_height + gap
        legend_y = 0

    if horizontal == 'right':
        new_width = max(main_width, legend_width)
        legend_x = new_width - legend_width
    else:  # left
        new_width = max(main_width, legend_width)
        legend_x = 0

    # Create the combined figure
    figure = Figure(new_width, new_height,
                    SVG(main_file).move(0, main_y),
                    SVG(legend_file).move(legend_x, legend_y))

    # Save the result
    figure.save(output_file)

def main():
    parser = argparse.ArgumentParser(description="Combine a main SVG with a legend SVG.")
    parser.add_argument("main_file", help="Path to the main SVG file")
    parser.add_argument("legend_file", help="Path to the legend SVG file")
    parser.add_argument("--vertical", choices=['top', 'bottom'], default='bottom',
                        help="Vertical position of the legend (default: bottom)")
    parser.add_argument("--horizontal", choices=['left', 'right'], default='right',
                        help="Horizontal position of the legend (default: right)")
    parser.add_argument("--output", default="combined_visualization.svg",
                        help="Output file name (default: combined_visualization.svg)")

    args = parser.parse_args()

    append_svg_legend(args.main_file, args.legend_file, args.vertical, args.horizontal, args.output)
    print(f"Combined SVG saved as {args.output}")

if __name__ == "__main__":
    main()
