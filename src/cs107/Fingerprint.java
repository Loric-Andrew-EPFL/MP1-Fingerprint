package cs107;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Provides tools to compare fingerprint.
 */
public class Fingerprint {

    /**
     * The number of pixels to consider in each direction when doing the linear
     * regression to compute the orientation.
     */
    public static final int ORIENTATION_DISTANCE = 16;

    /**
     * The maximum distance between two minutiae to be considered matching.
     */
    public static final int DISTANCE_THRESHOLD = 5;

    /**
     * The number of matching minutiae needed for two fingerprints to be considered
     * identical.
     */
    public static final int FOUND_THRESHOLD = 20;

    /**
     * The distance between two angle to be considered identical.
     */
    public static final int ORIENTATION_THRESHOLD = 20;

    /**
     * The offset in each direction for the rotation to test when doing the
     * matching.
     */
    public static final int MATCH_ANGLE_OFFSET = 2;

    /**
     * Creates a copy of a boolean 2D array
     *
     * @param image the array to be copied
     * @return an array with the same values and size as the argument but with different reference
     */
    public static boolean[][] createArrayCopy(boolean[][] image) {
        boolean[][] imageCopy = new boolean[image.length][image[0].length];
        for (int i = 0; i < image.length; i++) {
            for (int j = 0; j < image[0].length; j++) {
                imageCopy[i][j] = image[i][j];
            }
        }
        return imageCopy;
    }

    /**
     * Returns an array containing the value of the 8 neighbours of the pixel at
     * coordinates <code>(row, col)</code>.
     * <p>
     * The pixels are returned such that their indices corresponds to the following
     * diagram:<br>
     * ------------- <br>
     * | 7 | 0 | 1 | <br>
     * ------------- <br>
     * | 6 | _ | 2 | <br>
     * ------------- <br>
     * | 5 | 4 | 3 | <br>
     * ------------- <br>
     * <p>
     * If a neighbours is out of bounds of the image, it is considered white.
     * <p>
     * If the <code>row</code> or the <code>col</code> is out of bounds of the
     * image, the returned value should be <code>null</code>.
     *
     * @param image array containing each pixel's boolean value.
     * @param row   the row of the pixel of interest, must be between
     *              <code>0</code>(included) and
     *              <code>image.length</code>(excluded).
     * @param col   the column of the pixel of interest, must be between
     *              <code>0</code>(included) and
     *              <code>image[row].length</code>(excluded).
     * @return An array containing each neighbours' value.
     */
    public static boolean[] getNeighbours(boolean[][] image, int row, int col) {
        assert (image != null); // special case that is not expected (the image is supposed to have been checked earlier)

        boolean[] neighbours = new boolean[8];
        int[] rowDiff = {-1, -1, 0, 1, 1, 1, 0, -1};
        int[] colDiff = {0, 1, 1, 1, 0, -1, -1, -1};
        for (int i = 0; i < 8; ++i) {
            if (
                    row + rowDiff[i] >= 0 &&
                            row + rowDiff[i] < image.length &&
                            col + colDiff[i] >= 0 &&
                            col + colDiff[i] < image[0].length
            ) {
                neighbours[i] = image[row + rowDiff[i]][col + colDiff[i]];
            }
        }
        return neighbours;
    }

    /**
     * Computes the number of black (<code>true</code>) pixels among the neighbours
     * of a pixel.
     *
     * @param neighbours array containing each pixel value. The array must respect
     *                   the convention described in
     *                   {@link #getNeighbours(boolean[][], int, int)}.
     * @return the number of black neighbours.
     */
    public static int blackNeighbours(boolean[] neighbours) {
        int count = 0;
        for (boolean neighbour : neighbours) {
            if (neighbour) {
                count++;
            }
        }
        return count;
    }

    /**
     * Computes the number of white to black transitions among the neighbours of
     * pixel.
     *
     * @param neighbours array containing each pixel value. The array must respect
     *                   the convention described in
     *                   {@link #getNeighbours(boolean[][], int, int)}.
     * @return the number of white to black transitions.
     */
    public static int transitions(boolean[] neighbours) {
        int countTransitions = 0;
        for (int i = 0; i < neighbours.length; i++) {
            if (i == 0 && neighbours[0] && !neighbours[7]) {
                countTransitions++;
            } else if (i >= 1 && neighbours[i] && !neighbours[i - 1]) {
                countTransitions++;
            }
        }
        return countTransitions;
    }

    /**
     * Returns <code>true</code> if the images are identical and false otherwise.
     *
     * @param image1 array containing each pixel's boolean value.
     * @param image2 array containing each pixel's boolean value.
     * @return <code>True</code> if they are identical, <code>false</code>
     * otherwise.
     */
    public static boolean identical(boolean[][] image1, boolean[][] image2) {
        return Arrays.deepEquals(image1, image2);
    }

    /**
     * Internal method used by {@link #thin(boolean[][])}.
     *
     * @param image array containing each pixel's boolean value.
     * @param step  the step to apply, Step 0 or Step 1.
     * @return A new array containing each pixel's value after the step.
     */
    public static boolean[][] thinningStep(boolean[][] image, int step) {
        boolean[][] imageCopy = createArrayCopy(image);
        int[] positionIJ = step == 0 ? new int[]{4, 2} : new int[]{6, 0};
        for (int i = 0; i < image.length; i++) {
            for (int j = 0; j < image[0].length; j++) {
                boolean[] neighbours = getNeighbours(image, i, j);
                if (
                        image[i][j] &&
                                neighbours.length >= 1 &&
                                blackNeighbours(neighbours) >= 2 &&
                                blackNeighbours(neighbours) <= 6 &&
                                transitions(neighbours) == 1 &&
                                (!neighbours[0] /*P0*/ ||
                                        !neighbours[2] /*P2*/ ||
                                        !neighbours[positionIJ[0]] /*P4 or P6*/) &&
                                (!neighbours[positionIJ[1]] /*P2 or P0*/ ||
                                        !neighbours[4] /*P4*/ ||
                                        !neighbours[6] /*P6*/)
                ) {
                    imageCopy[i][j] = false;
                }
            }
        }
        return imageCopy;
    }

    /**
     * Compute the skeleton of a boolean image.
     *
     * @param image array containing each pixel's boolean value.
     * @return array containing the boolean value of each pixel of the image after
     * applying the thinning algorithm.
     */
    public static boolean[][] thin(boolean[][] image) {
        boolean[][] thinned = createArrayCopy(image);
        boolean[][] previousThinned;
        do {
            previousThinned = createArrayCopy(thinned);
            thinned = thinningStep(thinned, 0);
            thinned = thinningStep(thinned, 1);
        } while (!identical(thinned, previousThinned));

        return thinned;
    }

    /**
     * Computes all pixels that are connected to the pixel at coordinate
     * <code>(row, col)</code> and within the given distance of the pixel.
     *
     * @param image    array containing each pixel's boolean value.
     * @param row      the first coordinate of the pixel of interest.
     * @param col      the second coordinate of the pixel of interest.
     * @param distance the maximum distance at which a pixel is considered.
     * @return An array where <code>true</code> means that the pixel is within
     * <code>distance</code> and connected to the pixel at
     * <code>(row, col)</code>.
     */
    public static boolean[][] connectedPixels(boolean[][] image, int row, int col, int distance) {
        boolean[][] connectedPixels = new boolean[image.length][image[0].length];
        ArrayList<int[]> previousConnected = new ArrayList<>();

        int[] rowDiff = {-1, -1, 0, 1, 1, 1, 0, -1};
        int[] colDiff = {0, 1, 1, 1, 0, -1, -1, -1};

        int currentRow = row;
        int currentCol = col;

        previousConnected.add(new int[]{currentRow, currentCol});

        do {
            connectedPixels[currentRow][currentCol] = true;
            boolean[] currentNeighbours = getNeighbours(image, currentRow, currentCol);
            boolean[] connectedNeighbours = getNeighbours(connectedPixels, currentRow, currentCol);
            if (blackNeighbours(currentNeighbours) - blackNeighbours(connectedNeighbours) >= 1) {
                for (int i = 0; i < currentNeighbours.length; ++i) {
                    if (currentNeighbours[i] && !connectedNeighbours[i]) {
                        if (
                                Math.abs(currentRow + rowDiff[i] - row) <= distance &&
                                        Math.abs(currentCol + colDiff[i] - col) <= distance
                        ) {
                            currentRow += rowDiff[i];
                            currentCol += colDiff[i];
                            previousConnected.add(new int[]{currentRow, currentCol});
                            break;
                        } else if (i == currentNeighbours.length - 1) {
                            previousConnected.remove(previousConnected.size() - 1);
                            currentRow = previousConnected.get((previousConnected.size() - 1))[0];
                            currentCol = previousConnected.get((previousConnected.size() - 1))[1];
                        }
                    } else if (i == currentNeighbours.length - 1) {
                        previousConnected.remove(previousConnected.size() - 1);
                        currentRow = previousConnected.get((previousConnected.size() - 1))[0];
                        currentCol = previousConnected.get((previousConnected.size() - 1))[1];
                    }
                }
            } else {
                previousConnected.remove(previousConnected.size() - 1);
                if (previousConnected.size() > 0) {
                    currentRow = previousConnected.get((previousConnected.size() - 1))[0];
                    currentCol = previousConnected.get((previousConnected.size() - 1))[1];
                }
            }
        } while (previousConnected.size() > 0);

        return connectedPixels;
    }

    /**
     * Computes the slope of a minutia using linear regression.
     *
     * @param connectedPixels the result of
     *                        {@link #connectedPixels(boolean[][], int, int, int)}.
     * @param row             the row of the minutia.
     * @param col             the col of the minutia.
     * @return the slope.
     */
    public static double computeSlope(boolean[][] connectedPixels, int row, int col) {
        double sumOfXSquared = 0;
        double sumOfYSquared = 0;
        double sumOfXY = 0;

        for (int i = 0; i < connectedPixels.length; ++i) {
            for (int j = 0; j < connectedPixels[0].length; ++j) {
                if (connectedPixels[i][j]) {
                    sumOfXSquared += Math.pow(j - col, 2);
                    sumOfYSquared += Math.pow(row - i, 2);
                    sumOfXY += (j - col) * (row - i);
                }
            }
        }

        if (sumOfXSquared == 0) {
            return Double.POSITIVE_INFINITY;
        } else if (sumOfXSquared >= sumOfYSquared) {
            return sumOfXY / sumOfXSquared;
        } else {
            return sumOfYSquared / sumOfXY;
        }
    }

    /**
     * Computes the orientation of a minutia in radians.
     *
     * @param connectedPixels the result of
     *                        {@link #connectedPixels(boolean[][], int, int, int)}.
     * @param row             the row of the minutia.
     * @param col             the col of the minutia.
     * @param slope           the slope as returned by
     *                        {@link #computeSlope(boolean[][], int, int)}.
     * @return the orientation of the minutia in radians.
     */
    public static double computeAngle(boolean[][] connectedPixels, int row, int col, double slope) {
        int aboveCount = 0;
        int belowCount = 0;
        double angle = Math.atan(slope);

        for (int i = 0; i < connectedPixels.length; ++i) {
            for (int j = 0; j < connectedPixels[0].length; ++j) {
                if (connectedPixels[i][j] && (i != row && col != j)) {
                    double xPos = j - col;
                    double yPos = row - i;
                    if (slope != Double.POSITIVE_INFINITY) {
                        if (yPos >= ((-1) * xPos) / slope) {
                            aboveCount++;
                        } else {
                            belowCount++;
                        }
                    } else {
                        if (i > row) {
                            belowCount++;
                        } else {
                            aboveCount++;
                        }
                    }
                }
            }
        }
        if (slope != Double.POSITIVE_INFINITY) {
            if ((angle >= 0 && belowCount > aboveCount) || (angle <= 0 && aboveCount > belowCount)) {
                angle += Math.PI;
            }
        } else {
            if (belowCount > aboveCount) {
                angle *= (-1);
            }
        }

        return angle;
    }

    /**
     * Computes the orientation of the minutia that the coordinate <code>(row,
     * col)</code>.
     *
     * @param image    array containing each pixel's boolean value.
     * @param row      the first coordinate of the pixel of interest.
     * @param col      the second coordinate of the pixel of interest.
     * @param distance the distance to be considered in each direction to compute
     *                 the orientation.
     * @return The orientation in degrees.
     */
    public static int computeOrientation(boolean[][] image, int row, int col, int distance) {
        boolean[][] connectedPixels = connectedPixels(image, row, col, distance);
        double slope = computeSlope(connectedPixels, row, col);
        double angle = computeAngle(connectedPixels, row, col, slope);
        int angleDegrees = (int) Math.round((angle * 180) / Math.PI);
        while (angleDegrees < 0) {
            angleDegrees += 360;
        }
        while (angleDegrees > 359) {
            angleDegrees -= 360;
        }
        return angleDegrees;
    }

    /**
     * Extracts the minutiae from a thinned image.
     *
     * @param image array containing each pixel's boolean value.
     * @return The list of all minutiae. A minutia is represented by an array where
     * the first element is the row, the second is column, and the third is
     * the angle in degrees.
     * @see #thin(boolean[][])
     */
    public static List<int[]> extract(boolean[][] image) {
        List<int[]> minutiae = new ArrayList<>();
        for (int i = 1; i < image.length - 1; ++i) {
            for (int j = 1; j < image[0].length - 1; ++j) {
                if (image[i][j]) {
                    int transitions = transitions(getNeighbours(image, i, j));
                    if (transitions == 1 || transitions == 3) {
                        int orientation = computeOrientation(image, i, j, ORIENTATION_DISTANCE);
                        minutiae.add(new int[]{i, j, orientation});
                    }
                }
            }
        }
        return minutiae;
    }

    /**
     * Applies the specified rotation to the minutia.
     *
     * @param minutia   the original minutia.
     * @param centerRow the row of the center of rotation.
     * @param centerCol the col of the center of rotation.
     * @param rotation  the rotation in degrees.
     * @return the minutia rotated around the given center.
     */
    public static int[] applyRotation(int[] minutia, int centerRow, int centerCol, int rotation) {
        int[] rotatedMinutia = new int[3];
        int x = minutia[1] - centerCol;
        int y = centerRow - minutia[0];
        double cos = Math.cos((rotation * Math.PI) / 180);
        double sin = Math.sin((rotation * Math.PI) / 180);
        int newX = (int) Math.round(x * cos - y * sin);
        int newY = (int) Math.round(x * sin + y * cos);
        rotatedMinutia[0] = centerRow - newY;
        rotatedMinutia[1] = newX + centerCol;
        rotatedMinutia[2] = (minutia[2] + rotation) % 360;
        return rotatedMinutia;
    }

    /**
     * Applies the specified translation to the minutia.
     *
     * @param minutia        the original minutia.
     * @param rowTranslation the translation along the rows.
     * @param colTranslation the translation along the columns.
     * @return the translated minutia.
     */
    public static int[] applyTranslation(int[] minutia, int rowTranslation, int colTranslation) {
        int[] translatedMinutia = new int[3];
        translatedMinutia[0] = minutia[0] - rowTranslation;
        translatedMinutia[1] = minutia[1] - colTranslation;
        translatedMinutia[2] = minutia[2];
        return translatedMinutia;
    }

    /**
     * Computes the row, column, and angle after applying a transformation
     * (translation and rotation).
     *
     * @param minutia        the original minutia.
     * @param centerCol      the column around which the point is rotated.
     * @param centerRow      the row around which the point is rotated.
     * @param rowTranslation the vertical translation.
     * @param colTranslation the horizontal translation.
     * @param rotation       the rotation.
     * @return the transformed minutia.
     */
    public static int[] applyTransformation(int[] minutia, int centerRow, int centerCol, int rowTranslation,
                                            int colTranslation, int rotation) {
        int[] transformedMinutia = applyRotation(minutia, centerRow, centerCol, rotation);
        transformedMinutia = applyTranslation(transformedMinutia, rowTranslation, colTranslation);

        return transformedMinutia;
    }

    /**
     * Computes the row, column, and angle after applying a transformation
     * (translation and rotation) for each minutia in the given list.
     *
     * @param minutiae       the list of minutiae.
     * @param centerCol      the column around which the point is rotated.
     * @param centerRow      the row around which the point is rotated.
     * @param rowTranslation the vertical translation.
     * @param colTranslation the horizontal translation.
     * @param rotation       the rotation.
     * @return the list of transformed minutiae.
     */
    public static List<int[]> applyTransformation(List<int[]> minutiae, int centerRow, int centerCol, int rowTranslation,
                                                  int colTranslation, int rotation) {
        List<int[]> transformedMinutiae = new ArrayList<>();
        for (int[] minutia : minutiae) {
            transformedMinutiae.add(applyTransformation(minutia, centerRow, centerCol, rowTranslation, colTranslation,
                    rotation));
        }
        return transformedMinutiae;
    }

    /**
     * Counts the number of overlapping minutiae.
     *
     * @param minutiae1      the first set of minutiae.
     * @param minutiae2      the second set of minutiae.
     * @param maxDistance    the maximum distance between two minutiae to consider
     *                       them as overlapping.
     * @param maxOrientation the maximum difference of orientation between two
     *                       minutiae to consider them as overlapping.
     * @return the number of overlapping minutiae.
     */
    public static int matchingMinutiaeCount(List<int[]> minutiae1, List<int[]> minutiae2, int maxDistance,
                                            int maxOrientation) {
        int count = 0;
        for (int[] minutia1 : minutiae1) {
            for (int[] minutia2 : minutiae2) {
                double euclideanDistance = Math.sqrt(Math.pow(minutia1[0] - minutia2[0], 2) + Math.pow(minutia1[1] - minutia2[1], 2));
                int angleDifference = Math.abs(minutia1[2] - minutia2[2]);
                if (euclideanDistance <= maxDistance && angleDifference <= maxOrientation) {
                    count++;
                    //System.out.println(count + ": [" + minutia1[0] + ", " + minutia1[1] + ", " + minutia1[2] + "], " +"[" + minutia2[0] + ", " + minutia2[1] + ", " + minutia2[2] + "]");
                }
            }
        }
        return count;
    }

    /**
     * Compares the minutiae from two fingerprints.
     *
     * @param minutiae1 the list of minutiae of the first fingerprint.
     * @param minutiae2 the list of minutiae of the second fingerprint.
     * @return Returns <code>true</code> if they match and <code>false</code>
     * otherwise.
     */
    public static boolean match(List<int[]> minutiae1, List<int[]> minutiae2) {
        int maxFound = 0;
        for (int[] minutia1 : minutiae1) {
            for (int[] minutia2 : minutiae2) {
                int rotation = minutia1[2] - minutia2[2];
                for (int i = rotation - MATCH_ANGLE_OFFSET; i <= rotation + MATCH_ANGLE_OFFSET; ++i) {
                    List<int[]> transformedMinutia2 = applyTransformation(minutiae2, minutia2[0], minutia2[1], minutia2[0] - minutia1[0], minutia2[1] - minutia1[1], i);
                    int foundMinutia = matchingMinutiaeCount(minutiae1, transformedMinutia2, DISTANCE_THRESHOLD, ORIENTATION_THRESHOLD);
                    maxFound = Math.max(foundMinutia, maxFound);
                    if (foundMinutia >= FOUND_THRESHOLD) {
                        System.out.print("Match count: " + maxFound + ". ");
                        return true;
                    }
                }
            }
        }
        System.out.print("Match count: " + maxFound + ". ");
        return false;
    }
}
