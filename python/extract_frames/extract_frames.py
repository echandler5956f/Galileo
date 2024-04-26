#!.venv/bin/python3

# Extract keyframes from video: 
# A simple Python script using OpenCV and scikit-learn to extract keyframes 
# from all mp4 files in the current directory

import os
import cv2
import numpy as np
from sklearn.cluster import KMeans

def extract_frames_from_mp4(mp4_file, num_frames=10):
    # Create a folder for each video
    video_name = os.path.splitext(mp4_file)[0]
    os.makedirs(video_name, exist_ok=True)

    # Read the video
    cap = cv2.VideoCapture(mp4_file)
    total_frames = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))

    # Initialize variables
    frames = []
    frame_indices = []

    # Loop through frames
    for i in range(total_frames):
        ret, frame = cap.read()
        if ret:
            gray_frame = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
            frames.append(gray_frame.flatten())
            frame_indices.append(i)
        else:
            print(f"Error reading frame {i} from {mp4_file}")

    # Perform K-means clustering
    kmeans = KMeans(n_clusters=num_frames, n_init=num_frames, random_state=0).fit(frames)

    # Find the indices of the frames that are closest to the cluster centers
    selected_indices = [frame_indices[i] for i in np.argmin(kmeans.transform(frames), axis=0)]

    # Sort the frame indices in ascending order
    selected_indices = np.sort(selected_indices)

    # Extract and save the selected frames
    for i, idx in enumerate(selected_indices):
        cap.set(cv2.CAP_PROP_POS_FRAMES, idx)
        ret, frame = cap.read()
        if ret:
            # Add padding to make frames the same size
            padded_frame = cv2.copyMakeBorder(frame, 0, 0, 50, 50, cv2.BORDER_CONSTANT, value=(0, 0, 0))
            frame_path = os.path.join(video_name, f"frame_{i}.png")
            cv2.imwrite(frame_path, padded_frame)
        else:
            print(f"Error reading frame {idx} from {mp4_file}")

    cap.release()

if __name__ == "__main__":
    current_directory = os.getcwd()
    mp4_files = [file for file in os.listdir(current_directory) if file.lower().endswith(".mp4")]
    for mp4_file in mp4_files:
        print(f"Extracting frames from {mp4_file}")
        extract_frames_from_mp4(mp4_file)
        print(f"Extracted frames from {mp4_file}")