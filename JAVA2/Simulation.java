import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.JSlider;
import java.awt.Color;


public class Simulation extends Canvas implements Runnable {

	private boolean running;
	private DT mydt = new DT();
	private Embed myembed;
	private int XCENTER, YCENTER;
	private int XSIZE = 800, YSIZE = 800;
	private int speed = 200;
    private static final int BETA_MIN = 0;
    private static final int BETA_MAX = +100;
    private static final int BETA_INIT = +30;
    private static final int ALPHA_MIN = 0;
    private static final int ALPHA_MAX = +20;
    private static final int ALPHA_INIT = +5;

	public static void main(String[] args) {
		new Simulation();
		System.out.println("Launched simulation");
	}

	public Simulation() {
		XCENTER = XSIZE / 2;
		YCENTER = YSIZE / 2;
		setSize(XSIZE, YSIZE);
		Frame pictureFrame = new Frame("DT Simulation");
        pictureFrame.setBackground(Color.white);
		Panel canvasPanel = new Panel();
		canvasPanel.add(this);
		Panel controlPanel = new Panel();
		Button simControl = new Button("Start/Stop");

		simControl.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				running = !running;
			}
		});

        JLabel BETALabel = new JLabel("BETA*10", JLabel.CENTER);
        BETALabel.setAlignmentX(Component.CENTER_ALIGNMENT);
        JSlider sliderBETA = new JSlider(JSlider.HORIZONTAL,
                                        BETA_MIN, BETA_MAX, BETA_INIT);

        sliderBETA.addChangeListener(new ChangeListener(){
            public void stateChanged(ChangeEvent e) {
                JSlider source = (JSlider)e.getSource();
                if (!source.getValueIsAdjusting()) {
                    mydt.BETA = (double)source.getValue()*0.1; // multiply by 0.1 to get fractional BETA;
                    System.out.println("BETA is "+mydt.BETA);
                }
            }

        });

        sliderBETA.setMajorTickSpacing(20);
        sliderBETA.setMinorTickSpacing(1);
        sliderBETA.setPaintTicks(true);
        sliderBETA.setPaintLabels(true);

				/*
        JLabel KappaLabel = new JLabel("ALPHA*20", JLabel.CENTER);
        KappaLabel.setAlignmentX(Component.CENTER_ALIGNMENT);
        JSlider sliderKappa = new JSlider(JSlider.HORIZONTAL,
                                         ALPHA_MIN, ALPHA_MAX, ALPHA_INIT);

        sliderKappa.addChangeListener(new ChangeListener(){
            public void stateChanged(ChangeEvent e) {
                JSlider source = (JSlider)e.getSource();
                if (!source.getValueIsAdjusting()) {
                    mydt.ALPHA = (double)source.getValue()*0.05; // multiply by 0.1 to get fractional ALPHA/kappa_0??
                    System.out.println("ALPHA is "+mydt.ALPHA);
                }
            }

        });

        sliderKappa.setMajorTickSpacing(10);
        sliderKappa.setMinorTickSpacing(2);
        sliderKappa.setPaintTicks(true);
        sliderKappa.setPaintLabels(true);
				*/

        JLabel TimeLabel = new JLabel("Sweep time", JLabel.CENTER);
        TimeLabel.setAlignmentX(Component.CENTER_ALIGNMENT);
        JSlider sliderTime = new JSlider(JSlider.HORIZONTAL,
                                          50, 500, 100);

        sliderTime.addChangeListener(new ChangeListener(){
            public void stateChanged(ChangeEvent e) {
                JSlider source = (JSlider)e.getSource();
                if (!source.getValueIsAdjusting()) {
                    speed = source.getValue();
                    System.out.println("Inverse speed "+speed);
                }
            }

        });

        sliderTime.setMajorTickSpacing(100);
        sliderTime.setMinorTickSpacing(50);
        sliderTime.setPaintTicks(true);
        sliderTime.setPaintLabels(true);



        controlPanel.add(BETALabel);
        controlPanel.add(sliderBETA);
        //controlPanel.add(KappaLabel);
        //controlPanel.add(sliderKappa);
		controlPanel.add(simControl);
        controlPanel.add(TimeLabel);
		controlPanel.add(sliderTime);

		pictureFrame.add(canvasPanel);
		pictureFrame.add(controlPanel, BorderLayout.SOUTH);
		pictureFrame.addWindowListener(new WindowAdapter() {
			public void windowClosing(WindowEvent e) {
				System.exit(0);
			}
		});
		pictureFrame.pack();
		pictureFrame.setVisible(true);
		Thread myThread = new Thread(this);
		myThread.start();
	}

	public void run() {
		mydt.thermalize();
		while (true) {
			if (running) {
				for (int i = 0; i < mydt.VOL; i++) {
					mydt.trial_change();
				}
				// System.out.println("about to call embedding");
				mydt.tidy();
 				mydt.relabelnodes();
                mydt.laplacian();
                //System.out.println("Number of simplices "+DT.pointer_number);
				myembed = new Embed();
				myembed.ComputeEmbedding();
				repaint();
			}
			try {
				Thread.sleep(speed);
			} catch (InterruptedException e) {
			}
		}
	}

	public void paint(Graphics g) {
		int xPos1, yPos1, xPos2, yPos2, itmp, i, j;
		double SCALE = 300.0;
		// System.out.println("running");
		g.drawOval(XCENTER - (int) SCALE, YCENTER - (int) SCALE, (int) (2 * SCALE), (int) (2 * SCALE));
		for (i = 1; i < DT.node_number; i++)
			try {
				if (myembed.mark[i] != (DT.boundary_length+1)) {
					xPos2 = (int) (SCALE * myembed.x[i]) + XCENTER;
					yPos2 = (int) (SCALE * myembed.y[i]) + YCENTER;
                    for(j=DT.nstart[i];j<DT.nstart[i+1];j++){
						itmp = DT.ncol[j];
						xPos1 = (int) (SCALE * myembed.x[itmp]) + XCENTER;
						yPos1 = (int) (SCALE * myembed.y[itmp]) + YCENTER;

						g.drawLine(xPos2, yPos2, xPos1, yPos1);
					}
				}
			} catch (NullPointerException e) {}


	}

}
