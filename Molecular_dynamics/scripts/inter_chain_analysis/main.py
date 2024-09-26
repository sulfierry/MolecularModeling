# python main.py --input_dir vmd_frames --output_file detailed_interactions.tsv --prevalence_threshold 50 --top 30
import argparse
from interactions import ProcessInteractions, ViewInteractions


def main():
    parser = argparse.ArgumentParser(description='Process and visualize protein chain interactions.')
    parser.add_argument('--input_dir', type=str, default='vmd_frames', help='Directory containing input .dat files')
    parser.add_argument('--output_file', type=str, default='detailed_interactions.tsv', help='Output file name')
    parser.add_argument('--prevalence_threshold', type=float, default=50, help='Prevalence threshold for interactions')
    parser.add_argument('--top', type=int, default=30, choices=[10, 20, 30, 40, 50], help='Number of top interactions to visualize')
    args = parser.parse_args()

    # Process interactions
    process_interactions = ProcessInteractions(args.input_dir, args.output_file, args.prevalence_threshold)
    process_interactions.run()

    # View interactions
    view_interactions = ViewInteractions(args.output_file, args.top)
    view_interactions.run()

if __name__ == '__main__':
    main()
