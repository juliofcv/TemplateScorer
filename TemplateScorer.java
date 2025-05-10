package com.gd.bioutils.scoring;

import com.gd.bioutils.read.ReadISO19794;
import com.gd.bioutils.entity.Minutiae;
import java.awt.geom.Point2D;
import java.util.*;

public class TemplateScorer {
    private static final double DEFAULT_RESOLUTION = 500.0;
    private static final double MIN_RESOLUTION_THRESHOLD = 100.0;
    private static final double MAX_RESOLUTION_THRESHOLD = 10000.0;
    private static final double MIN_QUALITY_WEIGHT = 0.01;
    private static final double TRIM_RATIO = 0.1;
    private static final int HIST_BINS = 12;
    private static final double BETA = 0.5;
    private static final int MAX_MEDIAN_ITER = 10;
    private static final double MEDIAN_TOL = 1e-6;

    public static double scoreTemplate(ReadISO19794 template) {
        List<Minutiae> raw = template.getMinutiaeData();
        if (raw == null || raw.isEmpty()) return 0.0;

        List<Double> weights = calculateWeights(raw);
        List<Double> angles = extractAngles(raw);
        List<Point2D> points = extractPoints(raw);

        Point2D centroid = calculateRobustCentroid(points, weights);
        double[] resolutions = validateResolution(template.getXRes(), template.getYRes());

        List<Point2D> normalizedPoints = normalizePoints(points, centroid, resolutions, weights);
        double anglePCA = computePCAAlignment(normalizedPoints);
        List<Point2D> alignedPoints = rotatePoints(normalizedPoints, -anglePCA);

        List<Double> mahalanobisDistances = calculateMahalanobisDistances(alignedPoints, weights);
        double globalScore = calculateTrimmedMean(mahalanobisDistances);
        double quadrantScore = calculateQuadrantScore(alignedPoints);
        double orientationEntropy = calculateOrientationEntropy(angles);

        return globalScore * (1 + BETA * orientationEntropy) * quadrantScore;
    }

    private static List<Double> calculateWeights(List<Minutiae> raw) {
        List<Double> weights = new ArrayList<>();
        for (Minutiae m : raw) {
            double weight = Math.max(m.getQuality() / 100.0, MIN_QUALITY_WEIGHT);
            weights.add(weight);
        }
        return weights;
    }

    private static List<Double> extractAngles(List<Minutiae> raw) {
        List<Double> angles = new ArrayList<>();
        for (Minutiae m : raw) angles.add((double) m.getAngle());
        return angles;
    }

    private static List<Point2D> extractPoints(List<Minutiae> raw) {
        List<Point2D> points = new ArrayList<>();
        for (Minutiae m : raw) points.add(new Point2D.Double(m.getX(), m.getY()));
        return points;
    }

    private static Point2D calculateRobustCentroid(List<Point2D> points, List<Double> weights) {
        return geometricMedian(points, weights);
    }

    private static double[] validateResolution(double xRes, double yRes) {
        if (xRes < MIN_RESOLUTION_THRESHOLD || xRes > MAX_RESOLUTION_THRESHOLD) xRes = DEFAULT_RESOLUTION;
        if (yRes < MIN_RESOLUTION_THRESHOLD || yRes > MAX_RESOLUTION_THRESHOLD) yRes = DEFAULT_RESOLUTION;
        return new double[]{xRes, yRes};
    }

    private static List<Point2D> normalizePoints(List<Point2D> points, Point2D centroid, double[] resolutions, List<Double> weights) {
        List<Point2D> normalized = new ArrayList<>();
        double scale = calculateScaleFactor(points, centroid, resolutions, weights);
        for (Point2D p : points) {
            double x = (p.getX() - centroid.getX()) / resolutions[0] / scale;
            double y = (p.getY() - centroid.getY()) / resolutions[1] / scale;
            normalized.add(new Point2D.Double(x, y));
        }
        return normalized;
    }

    private static double calculateScaleFactor(List<Point2D> points, Point2D centroid, double[] resolutions, List<Double> weights) {
        double sum = 0, totalWeight = 0;
        for (int i = 0; i < points.size(); i++) {
            Point2D p = points.get(i);
            double distance = p.distance(centroid.getX(), centroid.getY());
            sum += (distance / ((resolutions[0] + resolutions[1]) / 2)) * weights.get(i);
            totalWeight += weights.get(i);
        }
        return sum / totalWeight;
    }

    private static double computePCAAlignment(List<Point2D> pts) {
        double sumXX=0, sumXY=0, sumYY=0;
        for (Point2D p: pts) {
            double x=p.getX(), y=p.getY();
            sumXX+=x*x; sumXY+=x*y; sumYY+=y*y;
        }
        int n=pts.size();
        double covXX=sumXX/n, covXY=sumXY/n, covYY=sumYY/n;
        return 0.5 * Math.atan2(2*covXY, covXX-covYY);
    }

    private static List<Point2D> rotatePoints(List<Point2D> pts, double angle) {
        double cos=Math.cos(angle), sin=Math.sin(angle);
        List<Point2D> out=new ArrayList<>(pts.size());
        for (Point2D p: pts) {
            out.add(new Point2D.Double(p.getX()*cos - p.getY()*sin, p.getX()*sin + p.getY()*cos));
        }
        return out;
    }

    /**
     * Mediana geométrica ponderada (Weiszfeld).
     */
    private static Point2D geometricMedian(List<Point2D> pts, List<Double> weights) {
        Point2D median = pts.get(0);
        for (int it = 0; it < MAX_MEDIAN_ITER; it++) {
            double wsum = 0, xsum = 0, ysum = 0;
            for (int i = 0; i < pts.size(); i++) {
                Point2D p = pts.get(i);
                double w = weights.get(i);
                double d = median.distance(p);
                double inv = (d < MEDIAN_TOL) ? 1 : w / d;
                xsum += p.getX() * inv;
                ysum += p.getY() * inv;
                wsum += inv;
            }
            Point2D next = new Point2D.Double(xsum / wsum, ysum / wsum);
            if (median.distance(next) < MEDIAN_TOL) break;
            median = next;
        }
        return median;
    }

    private static double[][] covariance(List<Point2D> pts, List<Double> weights) {
        double sumW = weights.stream().mapToDouble(Double::doubleValue).sum();
        Point2D mean = geometricMedian(pts, weights);
        double c00 = 0, c01 = 0, c11 = 0;
        for (int i = 0; i < pts.size(); i++) {
            Point2D p = pts.get(i);
            double w = weights.get(i);
            double dx = p.getX() - mean.getX(), dy = p.getY() - mean.getY();
            c00 += w * dx * dx;
            c01 += w * dx * dy;
            c11 += w * dy * dy;
        }
        return new double[][]{{c00 / sumW, c01 / sumW}, {c01 / sumW, c11 / sumW}};
    }

    private static List<Double> calculateMahalanobisDistances(List<Point2D> pts, List<Double> weights) {
        double[][] cov = covariance(pts, weights);
        double det = cov[0][0] * cov[1][1] - cov[0][1] * cov[1][0];
        double inv00 = cov[1][1] / det;
        double inv01 = -cov[0][1] / det;
        double inv11 = cov[0][0] / det;
        List<Double> dists = new ArrayList<>();
        for (Point2D p : pts) {
            double x = p.getX(), y = p.getY();
            double d2 = inv00 * x * x + 2 * inv01 * x * y + inv11 * y * y;
            dists.add(Math.sqrt(Math.max(d2, 0)));
        }
        Collections.sort(dists);
        return dists;
    }

    private static double calculateTrimmedMean(List<Double> dists) {
        int n = dists.size();
        int t = (int) (n * TRIM_RATIO);
        int from = Math.min(t, n - 1), to = Math.max(n - t, from + 1);
        double sum = 0;
        for (int i = from; i < to; i++) sum += dists.get(i);
        return sum / (to - from);
    }

    private static double calculateQuadrantScore(List<Point2D> pts) {
        int k = 3;
        double minX = pts.stream().mapToDouble(Point2D::getX).min().orElse(0);
        double maxX = pts.stream().mapToDouble(Point2D::getX).max().orElse(0);
        double minY = pts.stream().mapToDouble(Point2D::getY).min().orElse(0);
        double maxY = pts.stream().mapToDouble(Point2D::getY).max().orElse(0);
        double dx = (maxX - minX) / k, dy = (maxY - minY) / k;
        double prod = 1;
        int count = 0;
        for (int i = 0; i < k; i++) {
            for (int j = 0; j < k; j++) {
                double x0 = minX + i * dx, y0 = minY + j * dy;
                List<Double> block = new ArrayList<>();
                for (Point2D p : pts) {
                    if (p.getX() >= x0 && p.getX() < x0 + dx && p.getY() >= y0 && p.getY() < y0 + dy) {
                        block.add(p.distance(0, 0));
                    }
                }
                if (!block.isEmpty()) {
                    Collections.sort(block);
                    prod *= calculateTrimmedMean(block);
                    count++;
                }
            }
        }
        return Math.pow(prod, 1.0 / Math.max(count, 1));
    }

    private static double calculateOrientationEntropy(List<Double> angles) {
        double[] hist = new double[HIST_BINS];
        for (double a : angles) {
            int idx = (int) (a / 360.0 * HIST_BINS) % HIST_BINS;
            hist[idx]++;
        }
        int n = angles.size();
        double ent = 0;
        for (double h : hist) {
            if (h > 0) {
                double p = h / n;
                ent -= p * (Math.log(p) / Math.log(2));
            }
        }
        return ent / (Math.log(HIST_BINS) / Math.log(2));
    }
}
