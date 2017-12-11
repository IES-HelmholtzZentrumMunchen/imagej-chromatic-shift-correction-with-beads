import Objects3D.Counter3D;
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.measure.Calibration;
import ij.plugin.*;
import ij.plugin.filter.PlugInFilter;
import ij.process.AutoThresholder;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import ij.text.TextWindow;

import java.util.Vector;
import java.io.File;

/**
 * Chromatic Shift Correction With Beads
 *
 * A plugin for ImageJ 1.x to correct the chromatic aberration in fluorescent microscopy images using beads.
 * During the acquisition, take a few images with the same microscope settings but with beads instead of sample.
 * These beads will be used to define a chromatic aberration field and thereafter used to correct the chromatic shift on
 * the sample images.
 *
 * The beads detection is performed using Kozubek's method (see Kozubeck et al. (2000) 'An efficient algorithm for measurement
 * and correction of chromatic aberrations in fluorescent microscopy', Journal of Microscopy vol. 200, pt. 3, pp. 206-217).
 *
 * The field interpolation/extrapolation is performed using Sandwell's method (see Sandwell (1987) 'Biharmonic spline interpolation
 * of GEOS-3 and SEASAT altimeter data', Geophysical research letters, vol. 14, n. 2, pp. 139-142).
 *
 * @author Julien Pontabry
 */
public class ChromaticShiftCorrectionWithBeads implements PlugInFilter {
    /**
     * Title of the plugin to be displayed on top of dialogs.
     */
    public static String m_pluginTitle = "Chromatic shift correction";

    /**
     * Flag to know if a field or a calibration file have been loaded
     */
    private boolean m_loadField = false;

	/**
	 * @see ij.plugin.filter.PlugInFilter#setup(java.lang.String, ij.ImagePlus)
	 */
	@Override
	public int setup(String arg, ImagePlus imp) {
        //
        // GUI
        //

        // Create GUI
        final MyGenericDialog gui = new MyGenericDialog(m_pluginTitle);

        gui.addPathField("Calibration(s)", "", 30, true, true, false);
        gui.addPathField("Input image(s)", "", 30, true, true, false);
        gui.addPathField("Output chromatic field", "", 30, false, true, true);
        gui.addPathField("Output image(s)", "", 30, true, false, false);

        String[] choices = { "First", "Second" };
        gui.addRadioButtonGroup("Channel to be corrected", choices, 2, 30, "First");

        gui.addNumericField("Bead's diameter (in image's unit)", 0.2, 2);

        gui.showDialog();

        // Check GUI events
        if(gui.wasCanceled())
        {
            return DONE;
        }

        // Get paths
        String calibrationPath = gui.getNextString();
        String       inputPath = gui.getNextString();
        String      outputPath = gui.getNextString();
        String     outputIPath = gui.getNextString();
        String   channelChoice = gui.getNextRadioButton();

        int channelToCorrect = (channelChoice.compareTo(choices[0]) == 0 ? 1 : 2);

        double beadSize = gui.getNextNumber();


        //
        // Checking user inputs
        //

        // Check inputs
        if (calibrationPath.isEmpty() || inputPath.isEmpty()) {
            IJ.error("The calibration and input paths should not be empty !");
            return DONE;
        }

        if (Double.compare(beadSize, 0.0) <= 0) {
            IJ.error("The bead's size cannot be lesser than or equal to 0 !");
            return DONE;
        }

        // Check if folder are empty or not
        File calibrationFolder = new File(calibrationPath);
        File       inputFolder = new File(inputPath);

        File[] listOfCalibrationFiles = calibrationFolder.listFiles();
        File[]       listOfInputFiles = inputFolder.listFiles();

        if (calibrationFolder.isDirectory() && (listOfCalibrationFiles == null || listOfCalibrationFiles.length == 0)) {
            IJ.error("The calibration folder is empty !");
            return DONE;
        }

        if (inputFolder.isDirectory() && (listOfInputFiles == null || listOfInputFiles.length == 0)) {
            IJ.error("The input folder is empty !");
            return DONE;
        }


        //
        // Computing chromatic shift field
        //

        // Open calibration file(s)
        Vector< ImagePlus > calibrationImages = new Vector< ImagePlus >();

        if (calibrationFolder.isDirectory()) { // files in a folder
            if (listOfCalibrationFiles != null) {
                for (File file : listOfCalibrationFiles) {
                    if (file.isFile() && !file.getName().startsWith(".")) {
                        calibrationImages.add(IJ.openImage(file.getPath()));

                        if (calibrationImages.lastElement().getNChannels() != 2) {
                            IJ.error("This plugin supports only 2-channels calibration images !");
                            return DONE;
                        }
                    }
                }
            }
        }
        else { // !calibrationFolder.isDirectory() // one file
            calibrationImages.add(IJ.openImage(calibrationFolder.getPath()));

            if (calibrationImages.lastElement().getNChannels() != 2 && calibrationImages.lastElement().getNChannels() != 3) {
                IJ.error("This plugin supports either 2-channels calibration images or a chromatic field !");
                return DONE;
            }
        }

        // Check the number of files
        if (calibrationImages.size() == 0) {
            IJ.error("There is no calibration file !");
            return DONE;
        }

        // When one file, check if it is a field (32-bits) or a calibration
        if (calibrationImages.size() == 1 && calibrationImages.get(0).getBitDepth() == 32) {
            m_loadField = true;
        }

        // Define output and working sizes
        Vector< Integer > outputSize = new Vector< Integer >();
        outputSize.add(calibrationImages.get(0).getWidth());
        outputSize.add(calibrationImages.get(0).getHeight());
        outputSize.add(calibrationImages.get(0).getNSlices());

        Vector< Integer > workingSize = new Vector< Integer >();
        workingSize.addAll(outputSize);

        // Load the field (if any) or estimate it from calibration
        ImagePlus field;

        if (m_loadField) {
            field = calibrationImages.get(0);
        }
        else { // m_loadField
            // Convert bead's size from microns to pixels
            Calibration cal = calibrationImages.get(0).getCalibration();
                   beadSize = (beadSize / cal.pixelWidth);

            // Make isotropic voxels in calibration images (only for 3D)
            if (calibrationImages.get(0).getNSlices() > 1) {
                // Find the final sizes
                Calibration calib = calibrationImages.get(0).getCalibration();

                double scaleFactorX = calib.pixelWidth  / calib.pixelWidth;
                double scaleFactorY = calib.pixelHeight / calib.pixelWidth;
                double scaleFactorZ = calib.pixelDepth  / calib.pixelWidth;

                int  width = (int)(calibrationImages.get(0).getWidth()   * scaleFactorX);
                int height = (int)(calibrationImages.get(0).getHeight()  * scaleFactorY);
                int  depth = (int)(calibrationImages.get(0).getNSlices() * scaleFactorZ);

                // Rescale calibration images
                for (int i = 0; i < calibrationImages.size(); i++) {
                    ImagePlus  calibration = calibrationImages.get(i);
                    String calibrationName = calibration.getTitle();
                    IJ.run(calibration, "Scale...", "x=" + scaleFactorX + " y=" + scaleFactorY + " z=" + scaleFactorZ + " width=" + width + " height=" + height + " depth=" + depth + " interpolation=Bicubic average title=Rescaled");

                    IJ.selectWindow("Rescaled");
                    calibration = IJ.getImage();
                    calibration.hide();
                    calibration.setTitle(calibrationName);
                    calibrationImages.set(i, calibration);
                }

                // Modify the output size accordingly
                workingSize.set(0, width);
                workingSize.set(1, height);
                workingSize.set(2, depth);
            }

            // Beads detection
            Vector< Point >  firstListOfPoints = new Vector< Point >();
            Vector< Point > secondListOfPoints = new Vector< Point >();

            this.DetectMicroBeads(calibrationImages, firstListOfPoints, secondListOfPoints, beadSize, 0.3, 21);

            if (firstListOfPoints.size() == 0 | secondListOfPoints.size() == 0) {
                IJ.error("Unable to find micro-beads in calibration images !");
                return DONE;
            }

            // Compute chromatic deformation field
            if (channelToCorrect == 1) { // First channel
                field = this.EstimateChromaticFieldWithScatteredData(secondListOfPoints, firstListOfPoints, workingSize);
            }
            else { // channelToCorrect == 2 // Second channel
                field = this.EstimateChromaticFieldWithScatteredData(firstListOfPoints, secondListOfPoints, workingSize);
            }

            if (!outputPath.isEmpty()) {
                IJ.save(field, outputPath);
            }
        }

        // Memory clean
        calibrationImages = null;
        calibrationFolder = null;


        //
        // Application of the chromatic field to input images
        //

        // Open input file(s)
        Vector< ImagePlus > inputImages = new Vector< ImagePlus >();
        Vector< String >     inputNames = new Vector< String >();

        if (inputFolder.isDirectory()) { // files in a folder
            if (listOfInputFiles != null) {
                for (File file : listOfInputFiles) {
                    if (file.isFile() && !file.getName().startsWith(".")) {
                        inputImages.add(IJ.openImage(file.getPath()));
                        inputNames.add(file.getName());

                        if (inputImages.lastElement().getNChannels() != 2) {
                            IJ.error("This plugin supports only 2-channels images !");
                            return DONE;
                        }

                        if (!m_loadField && (inputImages.lastElement().getWidth() != outputSize.get(0) ||
                                inputImages.lastElement().getHeight() != outputSize.get(1) ||
                                inputImages.lastElement().getNSlices() != outputSize.get(2))) {
                            IJ.error("The size of image(s) is not the same as calibration(s) !");
                            return DONE;
                        }
                    }
                }
            }
        }
        else { // !folder.isDirectory()
            inputImages.add(IJ.openImage(inputFolder.getPath()));
            inputNames.add(inputFolder.getPath());

            if (inputImages.lastElement().getNChannels() != 2) {
                IJ.error("This plugin supports only 2-channels images !");
                return DONE;
            }

            if (!m_loadField && (inputImages.lastElement().getWidth() != outputSize.get(0) ||
                    inputImages.lastElement().getHeight() != outputSize.get(1) ||
                    inputImages.lastElement().getNSlices() != outputSize.get(2))) {
                IJ.error("The size of image(s) is not the same as calibration !");
                return DONE;
            }
        }

        // Check the number of files
        if (inputImages.size() == 0) {
            IJ.error("There is no input file !");
            return DONE;
        }

        if (m_loadField) {
            outputSize.set(0, inputImages.get(0).getWidth());
            outputSize.set(1, inputImages.get(0).getHeight());
            outputSize.set(2, inputImages.get(0).getNSlices());
        }

        // Make isotropic voxels in calibration images (only for 3D)
        if (inputImages.get(0).getNSlices() > 1) {
            double scaleFactorX = (double)workingSize.get(0) / (double)outputSize.get(0);
            double scaleFactorY = (double)workingSize.get(1) / (double)outputSize.get(1);
            double scaleFactorZ = (double)workingSize.get(2) / (double)outputSize.get(2);

            // Rescale calibration images
            for (int i = 0; i < inputImages.size(); i++) {
                ImagePlus  input = inputImages.get(i);
                String inputName = input.getTitle();
                IJ.run(input, "Scale...", "x=" + scaleFactorX + " y=" + scaleFactorY + " z=" + scaleFactorZ + " width=" + workingSize.get(0) + " height=" + workingSize.get(1) + " depth=" + workingSize.get(2) + " interpolation=Bicubic average title=Rescaled");

                IJ.selectWindow("Rescaled");
                input = IJ.getImage();
                input.hide();
                input.setTitle(inputName);
                inputImages.set(i, input);
            }
        }

        // Apply the correction to input files
        Vector< ImagePlus > outputImages = this.CorrectChromaticShift(inputImages, field, channelToCorrect);

        // Rescale input images to the output size (only for 3D)
        if (outputImages.get(0).getNSlices() > 1) {
            double scaleFactorX = (double)outputSize.get(0) / (double)workingSize.get(0);
            double scaleFactorY = (double)outputSize.get(1) / (double)workingSize.get(1);
            double scaleFactorZ = (double)outputSize.get(2) / (double)workingSize.get(2);

            // Rescale calibration images
            for (int i = 0; i < outputImages.size(); i++) {
                ImagePlus  output = outputImages.get(i);
                String outputName = output.getTitle();
                IJ.run(output, "Scale...", "x=" + scaleFactorX + " y=" + scaleFactorY + " z=" + scaleFactorZ + " width=" + outputSize.get(0) + " height=" + outputSize.get(1) + " depth=" + outputSize.get(2) + " interpolation=Bicubic average title=Rescaled");

                IJ.selectWindow("Rescaled");
                output = IJ.getImage();
                output.hide();
                output.setTitle(outputName);
                outputImages.set(i, output);
            }
        }

        // Save corrected images
        if (outputIPath.isEmpty()) {
            outputIPath = inputPath;
        }

        if (outputImages.size() > 1 && inputNames.size() > 1) {
            for (int fileIndex = 0; fileIndex < outputImages.size() && fileIndex < inputNames.size(); fileIndex++) {
                String fileName = inputNames.get(fileIndex);
                String[] s = fileName.split("\\.");
                IJ.save(outputImages.get(fileIndex), outputIPath + fileName.substring(0,fileName.length()-1-s[s.length-1].length()) + "_cs-corrected.tif");
            }
        }
        else { // outputImages.size() == 1
            String fileName = inputNames.get(0);
            String[]      s = fileName.split("\\.");
            IJ.save(outputImages.get(0), fileName.substring(0,fileName.length()-1-s[s.length-1].length()) + "_cs-corrected.tif");
        }

        // Memory clean
        inputImages  = null;
        outputImages = null;
        outputSize = null;

		return DONE;
	}

	/**
	 * @see ij.plugin.filter.PlugInFilter#run(ij.process.ImageProcessor)
	 */
	@Override
	public void run(ImageProcessor ip) {
        // ----
	}

    /**
     * Correct the chromatic shift with a calibration field.
     * @param images Input images.
     * @param field Input deformation field modeling chromatic shift.
     * @param channelToCorrect Channel to move with the deformation field.
     * @return Corrected images.
     */
    public Vector< ImagePlus > CorrectChromaticShift(Vector< ImagePlus > images, ImagePlus field, int channelToCorrect) {
        // Pre-condition
        assert                images != null;
        assert          images.size() > 0;
        assert                 field != null;
        assert      field.getWidth() == images.get(0).getWidth();
        assert     field.getHeight() == images.get(0).getHeight();
        assert    field.getNSlices() == images.get(0).getNSlices();
        assert (field.getNChannels() == 2 || field.getNChannels() == 3);
        assert  field.getNChannels() == images.get(0).getNChannels();

        // Initialise output
        Vector< ImagePlus > outputImages = new Vector< ImagePlus >();

        // Get field stack
        ImageStack fieldStack = field.getImageStack();

        // Correct each input image
        for (ImagePlus inputImage : images) {
            // Get input image stacks
            ImageStack inputStack = inputImage.getImageStack();

            // Create output image
            ImagePlus outputImage = inputImage.duplicate();

            // Get output image stacks
            ImageStack outputStack = outputImage.getImageStack();

            // Initialise output
            for (int z = 0; z < outputStack.getSize()/2; z++) {
                ImageProcessor proc = outputStack.getProcessor(z*2 + (channelToCorrect-1) + 1);
                proc.setValue(0);
                proc.fill();
            }

            // Go through all voxels
            for (int z = 0; z < outputStack.getSize()/2; z++) {
                for (int y = 0; y < outputStack.getHeight(); y++) {
                    for (int x = 0; x < outputStack.getWidth(); x++) {
                        outputStack.setVoxel(x,y,z*2+(channelToCorrect-1),
                                    this.linearlyInterpolate(inputStack, channelToCorrect-1,
                                            x + fieldStack.getVoxel(x,y,z*3),
                                            y + fieldStack.getVoxel(x,y,z*3+1),
                                            z + fieldStack.getVoxel(x,y,z*3+2)
                                    )
                                );
                    } // for x
                } // for y
            } // for z

            outputImages.add(outputImage);
        }

        // Post-conditions
        assert images.size() == outputImages.size();


        return outputImages;
    }

    /**
     * Get linearly interpolated value of continuous coordinates.
     * @param stack Image stack.
     * @param channelToCorrect Channel to correct in stack
     * @param x X-coordinate
     * @param y Y-coordinate
     * @param z Z-coordinate
     * @return Interpolated value.
     */
    private double linearlyInterpolate(ImageStack stack, int channelToCorrect, double x, double y, double z) {
        // Pre-compute sizes
        int  width = stack.getWidth();
        int height = stack.getHeight();
        int  depth = stack.getSize()/2;

        // Check for borders
        if (Double.compare(x,0) <= 0)
            x = 0;
        else if (Double.compare(x,width-1)  >= 0)
            x = width-1;

        if (Double.compare(y,0) <= 0)
            y = 0;
        else if (Double.compare(y,height-1)  >= 0)
            y = height-1;

        if (Double.compare(z,0) <= 0)
            z = 0;
        else if (Double.compare(z,depth-1)  >= 0)
            z = depth-1;

        // Compute surrounding coordinates
        double xBefore = Math.floor(x), xAfter = Math.floor(x + 1);
        double yBefore = Math.floor(y), yAfter = Math.floor(y + 1);
        double zBefore = Math.floor(z), zAfter = Math.floor(z + 1);

        if (Double.compare(x,width-1) >= 0) {
            xBefore = width - 2;
            xAfter  = width - 1;
        }

        if (Double.compare(y,height-1) >= 0) {
            yBefore = height - 2;
            yAfter  = height - 1;
        }

        if (Double.compare(z,depth-1) >= 0) {
            zBefore = depth - 2;
            zAfter  = depth - 1;
        }

        // Compute ratio position of asked point
        double rx = (x-xBefore) / (xAfter-xBefore);
        double ry = (y-yBefore) / (yAfter-yBefore);
        double rz = (z-zBefore) / (zAfter-zBefore);

        // Interpolate using trilinear method
        return this.trilinearInterpolation(rx, ry, rz,
                stack.getVoxel((int) xBefore, (int) yBefore, (int) (zBefore*2+channelToCorrect)),
                stack.getVoxel((int) xBefore, (int) yAfter,  (int) (zBefore*2+channelToCorrect)),
                stack.getVoxel((int) xAfter,  (int) yBefore, (int) (zBefore*2+channelToCorrect)),
                stack.getVoxel((int) xAfter,  (int) yAfter,  (int) (zBefore*2+channelToCorrect)),
                stack.getVoxel((int) xBefore, (int) yBefore, (int) (zAfter*2+channelToCorrect)),
                stack.getVoxel((int) xBefore, (int) yAfter,  (int) (zAfter*2+channelToCorrect)),
                stack.getVoxel((int) xAfter,  (int) yBefore, (int) (zAfter*2+channelToCorrect)),
                stack.getVoxel((int) xAfter,  (int) yAfter,  (int) (zAfter*2+channelToCorrect))
        );
    }

    /**
     * Trilinear interpolation.
     * @param rx Ratio position of the X-coordinate.
     * @param ry Ratio position of the Y-coordinate.
     * @param rz Ratio position of the Z-coordinate.
     * @param v_xm1_ym1_zm1 Value of the bottom left corner point (lower slice).
     * @param v_xm1_yp1_zm1 Value of the bottom right corner point (lower slice).
     * @param v_xp1_ym1_zm1 Value of the top left corner point (lower slice).
     * @param v_xp1_yp1_zm1 Value of the top right corner point (lower slice).
     * @param v_xm1_ym1_zp1 Value of the bottom left corner point (upper slice).
     * @param v_xm1_yp1_zp1 Value of the bottom right corner point (upper slice).
     * @param v_xp1_ym1_zp1 Value of the top left corner point (upper slice).
     * @param v_xp1_yp1_zp1 Value of the top right corner point (upper slice).
     * @return Interpolated value
     */
    private double trilinearInterpolation(double rx, double ry, double rz,
                                          double v_xm1_ym1_zm1, double v_xm1_yp1_zm1, double v_xp1_ym1_zm1, double v_xp1_yp1_zm1,
                                          double v_xm1_ym1_zp1, double v_xm1_yp1_zp1, double v_xp1_ym1_zp1, double v_xp1_yp1_zp1) {
        double C0 = this.bilinearInterpolation(rx, ry, v_xm1_ym1_zm1, v_xm1_yp1_zm1, v_xp1_ym1_zm1, v_xp1_yp1_zm1);
        double C1 = this.bilinearInterpolation(rx, ry, v_xm1_ym1_zp1, v_xm1_yp1_zp1, v_xp1_ym1_zp1, v_xp1_yp1_zp1);

        return this.linearInterpolation(rz, C0, C1);
    }

    /**
     * Bilinear interpolation.
     * @param rx Ratio position of the X-coordinate.
     * @param ry Ratio position of the Y-coordinate.
     * @param v_xm1_ym1 Value of the bottom left corner point.
     * @param v_xm1_yp1 Value of the bottom right corner point.
     * @param v_xp1_ym1 Value of the top left corner point.
     * @param v_xp1_yp1 Value of the top right corner point.
     * @return Interpolated value.
     */
    private double bilinearInterpolation(double rx, double ry,
                                         double v_xm1_ym1, double v_xm1_yp1, double v_xp1_ym1, double v_xp1_yp1) {
        double C0 = this.linearInterpolation(rx, v_xm1_ym1, v_xp1_ym1);
        double C1 = this.linearInterpolation(rx, v_xm1_yp1, v_xp1_yp1);

        return this.linearInterpolation(ry, C0, C1);
    }

    /**
     * Linear interpolation.
     * @param rx Ratio position of the X-coordinate.
     * @param v_xm1 Value of the previous point.
     * @param v_xp1 Value of the next point.
     * @return Interpolated value.
     */
    private double linearInterpolation(double rx, double v_xm1, double v_xp1) {
        return v_xm1 * (1-rx) + v_xp1 * rx;
    }

    /**
     * Detect the micro beads inside a set of calibration images.
     * @param calibrationImages A vector of calibration images.
     * @param firstListOfPoints A list of spatial locations in the first channel.
     * @param secondListOfPoints A list of spatial locations in the second channel.
     * @param beadSize Bead's size in pixels.
     * @param beadRadiusTolerance Tolerance when filtering with size in percent of beadRadius.
     * @param neighbourhoodSize Size of the search neighbourhood around beads (should be more than the expected shift).
     */
    public void DetectMicroBeads(Vector< ImagePlus > calibrationImages, Vector< Point > firstListOfPoints, Vector< Point > secondListOfPoints, double beadSize, double beadRadiusTolerance, int neighbourhoodSize) {
        for (ImagePlus calibrationImage : calibrationImages) {
            // Split channels
            ImagePlus[] channels = ChannelSplitter.split(calibrationImage);

            // Get the Rois centres from the first channel
            float[][] roisCentres = this.ExtractROIsFromImage(channels[0], beadSize, beadRadiusTolerance, neighbourhoodSize);

            // Get the micro beads from the first channel
            firstListOfPoints.addAll(this.ExtractBeadsCentresFromImage(channels[0], roisCentres, neighbourhoodSize));

            // Get the micro beads from the second channel
            secondListOfPoints.addAll(this.ExtractBeadsCentresFromImage(channels[1], roisCentres, neighbourhoodSize));
        }
    }

    /**
     * Extract the Rois in image around beads.
     * @param beadsImage Image of beads.
     * @param beadSize Bead's size in pixels.
     * @param beadRadiusTolerance Tolerance when filtering with size.
     * @param neighbourhoodSize Size of the search neighbourhood around beads (should be more than the expected shift).
     * @return A list of rois' centres.
     */
    private float[][] ExtractROIsFromImage(ImagePlus beadsImage, double beadSize, double beadRadiusTolerance,int neighbourhoodSize) {
        // Initialise neighbourhood's size
        int  halfNeighbourhoodSize = neighbourhoodSize/2;
        int halfNeighbourhoodSizeZ = 0;
        int    neighbourhoodPixels = neighbourhoodSize * neighbourhoodSize;

        // Add the pixels count for 3D if it is a stack (assuming isotropic voxels)
        if (beadsImage.getNSlices() > 1) {
            neighbourhoodPixels   *= neighbourhoodSize;
            halfNeighbourhoodSizeZ = halfNeighbourhoodSize;
        }

        // Get initial segmentation
        // Do a maximal projection on Z axis
        ZProjector projector = new ZProjector(beadsImage);
        projector.setMethod(ZProjector.MAX_METHOD);
        projector.doProjection();
        ImagePlus beadsImageProjected = projector.getProjection();

        // Get the values of automatic segmentation of the maximal projection
        ImageProcessor thresholdProcessor = beadsImageProjected.getProcessor();
        thresholdProcessor.setAutoThreshold(AutoThresholder.Method.Yen, true);

        // Apply threshold on all of the slices
        ImageStack beadsStack = beadsImage.duplicate().getImageStack();

        for (int z = 1; z <= beadsStack.getSize(); z++) {
            beadsStack.getProcessor(z).threshold((int)thresholdProcessor.getMinThreshold());
            beadsStack.getProcessor(z).resetThreshold();
        }

        beadsImageProjected.close();

        // Extract centre of blobs
        double minBeadRadius = beadSize*(1.0-beadRadiusTolerance)/2.0;
        double maxBeadRadius = beadSize*(1.0+beadRadiusTolerance)/2.0;

        int minBeadPixels, maxBeadPixels;

        if (beadsImage.getNSlices() > 1) {
            minBeadPixels = (int) Math.round(4.0 * Math.PI * minBeadRadius * minBeadRadius * minBeadRadius / 3.0);
            maxBeadPixels = (int) Math.round(4.0 * Math.PI * maxBeadRadius * maxBeadRadius * maxBeadRadius / 3.0);
        }
        else { // beadsImage.getNSlices() == 1
            minBeadPixels = (int) Math.round(Math.PI * minBeadRadius * minBeadRadius);
            maxBeadPixels = (int) Math.round(Math.PI * maxBeadRadius * maxBeadRadius);
        }

        Counter3D counter = new Counter3D(new ImagePlus("TMP",beadsStack), 128, minBeadPixels, maxBeadPixels, true, false);
        float[][] blobsCentres = counter.getCentroidList();

        // Define Rois around beads' centres
        ImagePlus       roisImage = beadsImage.duplicate();
        ImageStack roisImageStack = roisImage.getImageStack();

        roisImage.deleteRoi();

        for (int z = 1; z <= roisImageStack.getSize(); z++) {
            ImageProcessor proc = roisImageStack.getProcessor(z);
            proc.setValue(0);
            proc.fill();
        }

        for (int centreIndex = 0; centreIndex < blobsCentres.length; centreIndex++) {
            int      x = Math.round(blobsCentres[centreIndex][0])-halfNeighbourhoodSize;
            int      y = Math.round(blobsCentres[centreIndex][1])-halfNeighbourhoodSize;
            int zBegin = Math.round(blobsCentres[centreIndex][2])-halfNeighbourhoodSizeZ;
            int   zEnd = Math.round(blobsCentres[centreIndex][2])+halfNeighbourhoodSizeZ;

            for (int z = zBegin; z <= zEnd; z++) { // For each slice
                roisImage.setPositionWithoutUpdate(1, z, 1);
                roisImage.setRoi(x, y, neighbourhoodSize, neighbourhoodSize);

                ImageProcessor proc = roisImageStack.getProcessor(z);
                proc.setValue(255);
                proc.fill(roisImage.getRoi());
            }
        }

        roisImage.deleteRoi();

        // Filter rois (keep only non-overlapping and non-touching rois)
        counter = new Counter3D(roisImage, 128, neighbourhoodPixels, neighbourhoodPixels, true, false);
        float[][] roisCentres = counter.getCentroidList();

        roisImage.close();

        // Correct the Z-coordinate for output rois centers
        for (int i = 0; i < roisCentres.length; i++) {
            roisCentres[i][2] -= 1;
        }

        return roisCentres;
    }

    /**
     * Extract beads' centres given a list of roi centres.
     * @param beadsImage Image of beads.
     * @param roisCentres List of rois' centres.
     * @param neighbourhoodSize Size of the search neighbourhood around beads (should be more than the expected shift).
     * @return Estimated centres of beads for the given image.
     */
    private Vector< Point > ExtractBeadsCentresFromImage(ImagePlus beadsImage, float[][] roisCentres, int neighbourhoodSize) {
        // Pre-conditions
        assert roisCentres.length > 0;
        assert beadsImage != null;

        // Initialise output
        Vector< Point > listOfPoints = new Vector< Point >();

        // Initialise neighbourhood's size
        int  halfNeighbourhoodSize = neighbourhoodSize/2;
        int halfNeighbourhoodSizeZ = 0;

        if (beadsImage.getNSlices() > 1) {
            halfNeighbourhoodSizeZ = halfNeighbourhoodSize;
        }

        // Estimate beads center positions
        for (int centreIndex = 0; centreIndex < roisCentres.length; centreIndex++) {
            // Pre-compute values
            int roi_first_y = Math.round(roisCentres[centreIndex][1]) - halfNeighbourhoodSize;
            int roi_first_x = Math.round(roisCentres[centreIndex][0]) - halfNeighbourhoodSize;
            int roi_begin_z = Math.round(roisCentres[centreIndex][2]) - halfNeighbourhoodSizeZ;
            int   roi_end_z = Math.round(roisCentres[centreIndex][2]) + halfNeighbourhoodSizeZ;

            // Get roi
            beadsImage.setRoi(roi_first_x, roi_first_y, neighbourhoodSize, neighbourhoodSize);
            Duplicator dup = new Duplicator();
            ImageStack currentRoi = dup.run(beadsImage, roi_begin_z+1, roi_end_z+1).getImageStack();

            // Get statistics from image intensities
            double medianStat = 0, maxStat = 0;

            for (int z = 1; z <= currentRoi.getSize(); z++) {
                ImageProcessor   proc = currentRoi.getProcessor(z);
                ImageStatistics stats = proc.getStatistics();

                if (stats.max > maxStat)
                    maxStat = stats.max;

                medianStat += stats.median;
            }

            medianStat /= currentRoi.getSize();

            // Compute a threshold the image using Kozubeck's formula
            int threshold = (int)Math.round((medianStat+maxStat)/2.0);

            // Go through each pixel to compute the center of bead
            Point beadCenter = new Point();
            int intensitiesAboveThreshold = 0;

            for (int z = 0; z < currentRoi.getSize(); z++) {
                for (int y = 0; y < currentRoi.getHeight(); y++) {
                    for (int x = 0; x < currentRoi.getWidth(); x++) {
                        double pixelValue = currentRoi.getVoxel(x,y,z);

                        if (pixelValue >= threshold) {
                            beadCenter.x += pixelValue * (x + roi_first_x);
                            beadCenter.y += pixelValue * (y + roi_first_y);
                            beadCenter.z += pixelValue * (z + roi_begin_z);

                            intensitiesAboveThreshold += pixelValue;
                        }
                    }
                }
            }

            beadCenter.x /= intensitiesAboveThreshold;
            beadCenter.y /= intensitiesAboveThreshold;
            beadCenter.z /= intensitiesAboveThreshold;

            listOfPoints.add(beadCenter);
        }

        return listOfPoints;
    }

    /**
     * Estimate the chromatic deformation field with two lists of points.
     * @param firstListOfPoints First list of points.
     * @param secondListOfPoints Second list of points.
     * @param outputSize Size of the output image (deformation field).
     * @return An image of the deformation field.
     */
    public ImagePlus EstimateChromaticFieldWithScatteredData(Vector< Point > firstListOfPoints, Vector< Point > secondListOfPoints, Vector< Integer > outputSize) {
        // Pre-conditions
        assert firstListOfPoints.size() == secondListOfPoints.size();

        // Compute the sparse deformation field
        Vector< Vector< Double > > listOfDeformationVectors = this.ComputeSparseDeformationField(firstListOfPoints, secondListOfPoints);

        // Remove too close points
        this.MergeClosePoints(firstListOfPoints, listOfDeformationVectors);
//TextWindow log = new TextWindow("pos", "", 500, 250);
//for (Point p : firstListOfPoints) {
//    log.append(p+"");
//}
//TextWindow log2 = new TextWindow("vec", "", 500, 250);
//for (Vector< Double > component : listOfDeformationVectors) {
//    String s = "";
//    for (Double value : component) {
//        s += value + " ";
//    }
//    log2.append(s);
//}
        // Reconstruct the whole field and return it
        return this.InterpolateScatteredDataWithBiharmonicSpline(firstListOfPoints, listOfDeformationVectors, outputSize);
    }

    /**
     * Compute a sparse deformation field with two lists of points.
     * @param firstPositions First list of points.
     * @param secondPositions Second list of points.
     * @return A sparse deformation field.
     */
    private Vector< Vector< Double > > ComputeSparseDeformationField(Vector< Point > firstPositions, Vector< Point > secondPositions) {
        // Pre-conditions
        assert firstPositions != null;
        assert secondPositions != null;
        assert firstPositions.size() == secondPositions.size();

        // Initialise the output
        Vector< Vector< Double > > sparseField = new Vector< Vector< Double > >();

        // Compute the displacement, component by component
        for (int componentIndex = 0; componentIndex < 3; componentIndex++) {
            sparseField.add(new Vector< Double >());

            for (int pointIndex = 0; pointIndex < firstPositions.size(); pointIndex++) {
                if (componentIndex == 0) {
                    sparseField.get(componentIndex).add(secondPositions.get(pointIndex).x - firstPositions.get(pointIndex).x);
                }
                else if (componentIndex == 1) {
                    sparseField.get(componentIndex).add(secondPositions.get(pointIndex).y - firstPositions.get(pointIndex).y);
                }
                else { // componentIndex != 0 && componentIndex != 1
                    sparseField.get(componentIndex).add(secondPositions.get(pointIndex).z - firstPositions.get(pointIndex).z);
                }
            } // for each point
        } // for each component

        // Post-conditions
        assert sparseField.size() == 2;
        assert sparseField.get(0).size() == firstPositions.size();

        return sparseField;
    }

    /**
     * Merge points that are close to each other (depending on the machine precision).
     * The positions and the values of points that are merged averaged.
     * @param positions Spatial positions of the points.
     * @param response Response of each point.
     */
    private void MergeClosePoints(Vector< Point > positions, Vector< Vector< Double > > response) {
        // Pre-condition
        assert positions != null;
        assert response != null;
        assert response.size() == 2;
        assert response.get(0) != null;
        assert positions.size() == response.get(0).size();

        // Find max and min for x,y and z components
        double minX = 0.0, maxX = 0.0, minY = 0.0, maxY = 0.0, minZ = 0.0, maxZ = 0.0;

        for (Point point : positions) {
            if (point.x < minX) {
                minX = point.x;
            }

            if (point.y < minY) {
                minY = point.y;
            }

            if (point.z < minZ) {
                minZ = point.z;
            }

            if (point.x > maxX) {
                maxX = point.x;
            }

            if (point.y > maxY) {
                maxY = point.y;
            }

            if (point.z > maxZ) {
                maxZ = point.z;
            }
        }

        // Pre-compute
        double epsilonx  = Math.pow(Math.ulp(0.5 * (maxX-minX)), 1.0/3.0);
        double epsilony  = Math.pow(Math.ulp(0.5 * (maxY-minY)), 1.0/3.0);
        double epsilonz  = Math.pow(Math.ulp(0.5 * (maxZ-minZ)), 1.0/3.0);
        double epsilonx2 = epsilonx * epsilonx;
        double epsilony2 = epsilony * epsilony;
        double epsilonz2 = epsilonz * epsilonz;
        double epsilon   = epsilonx2 + epsilony2 + epsilonz2;

        // Search for close points and remove them
        for (int firstPointIndex = 0; firstPointIndex < positions.size(); firstPointIndex++) {
            for (int secondPointIndex = 0; secondPointIndex < positions.size(); secondPointIndex++) {
                if (firstPointIndex != secondPointIndex) {
                    double diffx = positions.get(firstPointIndex).x - positions.get(secondPointIndex).x;
                    double diffy = positions.get(firstPointIndex).y - positions.get(secondPointIndex).y;
                    double diffz = positions.get(firstPointIndex).z - positions.get(secondPointIndex).z;

                    // If they are too close, average them
                    if (Double.compare(diffx * diffx + diffy * diffy + diffz * diffz, epsilon) <= 0) {
                        positions.get(firstPointIndex).x = 0.5 * (positions.get(firstPointIndex).x + positions.get(secondPointIndex).x);
                        positions.get(firstPointIndex).y = 0.5 * (positions.get(firstPointIndex).y + positions.get(secondPointIndex).y);
                        positions.get(firstPointIndex).z = 0.5 * (positions.get(firstPointIndex).z + positions.get(secondPointIndex).z);

                        response.get(0).set(firstPointIndex, 0.5 * (response.get(0).get(firstPointIndex) + response.get(0).get(secondPointIndex)));
                        response.get(1).set(firstPointIndex, 0.5 * (response.get(1).get(firstPointIndex) + response.get(1).get(secondPointIndex)));
                        response.get(2).set(firstPointIndex, 0.5 * (response.get(2).get(firstPointIndex) + response.get(2).get(secondPointIndex)));

                        positions.removeElementAt(secondPointIndex);
                        response.get(0).removeElementAt(secondPointIndex);
                        response.get(1).removeElementAt(secondPointIndex);
                        response.get(2).removeElementAt(secondPointIndex);
                    }
                } // else firstPointIndex == secondPointIndex
            } // for each second point
        } // for each first point
    }

    /**
     * Interpolate scattered data using a biharmonic spline method.
     * The method is from Sandwell (1987) 'Biharmonic spline interpolation
     * of GEOS-3 and SEASAT altimeter data', Geophysical research letters, vol. 14, n. 2, pp. 139-142.
     * @param positions Spatial locations of the scattered data.
     * @param response Response values corresponding the the spatial locations.
     * @param outputSize Size of image output.
     * @return An interpolated field over the output image.
     */
    private ImagePlus InterpolateScatteredDataWithBiharmonicSpline(Vector< Point > positions, Vector< Vector< Double > > response, Vector< Integer > outputSize) {
        // Pre-condition
        assert response != null;
        assert positions != null;
        assert response.get(0) != null;
        assert positions.size() == response.get(0).size() : "The number of points and corresponding response values should be the same.";
        assert outputSize.size() == 2 : "The number of image size should be 2 (2D images).";

        // TEST voir si c'est mieux avec la fonction de la bonne dimension
        // Compute the Green's function for data points (|p_i-p_j|^2 * (log(|p_i-p_j|)-1))
        double[][] greenFunction = new double[positions.size()][positions.size()];

        for (int i = 0; i < greenFunction.length; i++) {
            for (int j = 0; j < greenFunction[i].length; j++) {
                if (i != j) {
                    double    diffx = positions.get(i).x - positions.get(j).x;
                    double    diffy = positions.get(i).y - positions.get(j).y;
                    double    diffz = positions.get(i).z - positions.get(j).z;
                    double distance = Math.sqrt(diffx * diffx + diffy * diffy + diffz * diffz);

                    greenFunction[i][j] = (distance * distance) * (Math.log(distance) - 1.0);
                }
                else { // i == j
                    greenFunction[i][j] = 0.0;
                }
            } // for j
        } // for i

        // Initialise output
        Vector< ImagePlus > field = new Vector< ImagePlus >();

        try {
            // Solve the following linear system greenFunction.weights = response
            // to find the weights of the biharmonic spline function
            Vector< double[] > weights = new Vector< double[] >();

            for (int c = 0; c < response.size(); c++) {
                weights.add(this.SolveLinearSystem(greenFunction, response.get(c)));

                // Invariant
                assert weights.get(c) != null;
                assert weights.get(c).length == positions.size() : "The number of weights should be equal to the number of input points.";
            }

            // Initialise components of deformation fields
            for (int c = 0; c < weights.size(); c++) {
                field.add(IJ.createImage("Chromatic field C" + (c + 1), "32-bit black", outputSize.get(0), outputSize.get(1), 1, outputSize.get(2), 1));
            }

            Vector< ImageStack > components = new Vector< ImageStack >();

            for (int c = 0; c < field.size(); c++) {
                components.add(field.get(c).getImageStack());
            }

            // Pre-compute
            int  width = field.get(0).getWidth();
            int height = field.get(0).getHeight();
            int  depth = field.get(0).getNSlices();

            // Display a progress bar
            int numberOfVisitedVoxels = 0;
            int        numberOfVoxels =  width * height * depth;
            IJ.showProgress(numberOfVisitedVoxels, numberOfVoxels);
            IJ.showStatus("Interpolation of scattered deformation field...");

            // Go through all voxels
            for (int z = 0; z < depth; z++) {
                for (int y = 0; y < height; y++) {
                    for (int x = 0; x < width; x++) {
                        // Compute the Green's function for current position
                        double[] localGreenFunction = new double[positions.size()];

                        for (int i = 0; i < localGreenFunction.length; i++) {
                            double diffx = x - positions.get(i).x;
                            double diffy = y - positions.get(i).y;
                            double diffz = z - positions.get(i).z;
                            double distance = Math.sqrt(diffx * diffx + diffy * diffy + diffz * diffz);

                            localGreenFunction[i] = (distance * distance) * (Math.log(distance) - 1.0);
                        }

                        // Interpolate at current position
                        for (int c = 0; c < weights.size(); c++) {
                            double interpolationValue = 0.0;

                            for (int i = 0; i < localGreenFunction.length; i++) {
                                interpolationValue += localGreenFunction[i] * weights.get(c)[i];
                            }

                            components.get(c).setVoxel(x,y,z, (float) interpolationValue);
                        }

                        IJ.showProgress(++numberOfVisitedVoxels, numberOfVoxels);
                    }
                }
            }

            IJ.showProgress(numberOfVoxels, numberOfVoxels);
        } catch (Exception e) {
            IJ.showMessage("Error: " + e.getMessage());
        }

        // Merge components in a final deformation field
        RGBStackMerge merger = new RGBStackMerge();
        ImagePlus[]      tmp = new ImagePlus[field.size()];

        return merger.mergeHyperstacks(field.toArray(tmp), false);
    }

    /**
     * Solve a linear system by Gauss-Jordan triangulation.
     * The linear system to solve is Ax = b.
     * This method does not change the values of the inputs A and b.
     * @param A Coefficient matrix of the linear system.
     * @param b Response of the linear system.
     * @return The computed solution of the linear system.
     * @throws Exception If the matrix A is singular or near singular (or badly conditioned).
     */
    private double[] SolveLinearSystem(double[][] A, Vector< Double > b) throws Exception {
        // Pre-conditions
        assert A != null;
        assert A[0] != null;
        assert b != null;
        assert A.length == b.size() : "The number of rows of the matrix of equations should be equal to the size of the vector b.";

        int    numberOfRows = A.length;
        int numberOfColumns = A[0].length;

        // Make a deep copy of the A matrix (to prevent changes in this matrix)
        double[][] Acopy = new double[A.length][A[0].length];

        for (int i = 0; i < Acopy.length; i++) {
            for (int j = 0; j < Acopy[0].length; j++) {
                Acopy[i][j] = A[i][j];
            }
        }

        // First triangulate the matrix using the Gauss-Jordan algorithm
        for (int k = 0; k < numberOfRows; k++) {
            // Search for maximal pivot in column k
            double pivotValue = Acopy[k][k];
            int    pivotIndex = k;

            for (int i = k+1; i < numberOfRows; i++) {
                if (Math.abs(Acopy[i][k]) > pivotValue) {
                    pivotValue = Acopy[i][k];
                    pivotIndex = i;
                }
            }

            if (Double.compare(pivotValue,0.0) == 0) {
                throw new Exception("Matrix is singular or near singular !");
            }

            // Permute lines if it is needed
            if (pivotIndex != k) {
                for (int j = k; j < numberOfColumns; j++) {
                    double tmpValue = Acopy[k][j];
                    Acopy[k][j] = Acopy[pivotIndex][j];
                    Acopy[pivotIndex][j] = tmpValue;
                }

                double tmpValue = b.get(k);
                b.set(k, b.get(pivotIndex));
                b.set(pivotIndex, tmpValue);
            }

            // Apply the Gauss-Jordan reduction with pivot k
            for (int i = k+1; i < numberOfRows; i++) {
                double coefficient = Acopy[i][k] / pivotValue;

                for (int j = k; j < numberOfColumns; j++) {
                    Acopy[i][j] -= Acopy[k][j] * coefficient;
                }

                b.set(i, b.get(i) - b.get(k)*coefficient);
            }
        }

        // Initialise the solution vector
        double[] x = new double[numberOfRows];

        // Then compute the solution of this superior triangular matrix
        x[numberOfRows - 1] = b.get(numberOfRows-1) / Acopy[numberOfRows-1][numberOfRows-1];

        for (int i = numberOfRows-2; i >= 0; i--) {
            double sum = 0.0;

            for (int k = i+1; k < numberOfRows; k++) {
                sum += Acopy[i][k]*x[k];
            }

            x[i] = (b.get(i)-sum)/Acopy[i][i];
        }

        // Post-condition
        assert numberOfRows == x.length : "The number of rows of the matrix of equations and the size of the solution vector should be equal.";


        return x;
    }

	/**
	 * Main method for debugging.
	 *
	 * For debugging, it is convenient to have a method that starts ImageJ, loads an
	 * image and calls the plugin, e.g. after setting breakpoints.
	 *
	 * @param args unused
	 */
	public static void main(String[] args) {
		// set the plugins.dir property to make the plugin appear in the Plugins menu
		Class<?> clazz = ChromaticShiftCorrectionWithBeads.class;
		String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
		String pluginsDir = url.substring(5, url.length() - clazz.getName().length() - 6);
		System.setProperty("plugins.dir", pluginsDir);

		// start ImageJ
		new ImageJ();

        // run the plugin
        IJ.runPlugIn(clazz.getName(), "");
	}
}

class Point {
    /**
     * Coordinates of the 3D point (for 2D, z=0).
     */
    public double x,y,z;

    /**
     * Default constructor.
     * Initialise the point to null.
     */
    public Point() {
        x = 0.0;
        y = 0.0;
        z = 0.0;
    }

    /**
     * Constructor for 2D.
     * @param _x X coordinate.
     * @param _y Y coordinate.
     */
    public Point(double _x, double _y) {
        x = _x;
        y = _y;
    }

    /**
     * Constructor for 3D.
     * @param _x X coordinate.
     * @param _y Y coordinate.
     * @param _z Z coordinate.
     */
    public Point(double _x, double _y, double _z) {
        x = _x;
        y = _y;
        z = _z;
    }

    /**
     * Copy constructor.
     * @param p point.
     */
    public Point(Point p) {
        x = p.x;
        y = p.y;
        z = p.z;
    }

    public String toString() {
        return "(" + x + "," + y + "," + z + ")";
    }
}
