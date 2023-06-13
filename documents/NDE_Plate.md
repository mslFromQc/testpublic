NDE Open File Format

NDE 3.1.0 File Format: Plate Conventions

![Logo Description automatically
generated](C:/Dev/nde-format/documents/media/image1.png){width="1.885in"
height="2.5in"}

Evident

Version 0.9

April 1, 2023

# NDE File Format: Plate Conventions 

-   By the .nde file format's convention, surface **id:0** is the outer
    surface (or top for plates) and surface **id:1** is the inner
    surface (or bottom for plates).

-   **Units:** Angles are expressed in degrees, while all other units in
    the dataset are in International System (SI) units, expressed in
    meters and seconds unless otherwise indicated by a **\"unit\"** key.
    For example, the \"Bitfield\" and \"Percent\" units are used for the
    **ascan** in the dataset object.

## Axis and Coordinate System Definitions

-   **U**, **V**, **W**: Local or **surface** coordinate system. It is
    the main coordinate by which everything is related in each file. It
    describes positions on the specimen surface (**u**, **v**) and depth
    (**w**).

    -   **origin**: Set on the surface of the specimen and some
        preferred direction along the surface. For example, **U** can be
        set along the center line of a weld.

    -   The **U** axis is defined by the scenario. In a general mapping
        scenario, it can be described as the **\"active\"** axis, which
        often corresponds to the **Scan** axis and is the one on which
        probes are usually swept in one-line scan acquisition. For a
        weld scenario, **U** is defined along the weld.

    -   The **V** axis is perpendicular to **U**. In most cases, it
        corresponds to the **Index** axis in raster acquisition.

    -   The **W** axis is normal to the surface and points inside the
        material. **W** is particularly useful in total focusing method
        (TFM) acquisition to position the data. However, in most cases,
        an **\"Ultrasound\"** axis is used to describe the data instead.

    -   The **U** and **V** axis properties are given in the
        **dataEncodings** object of the domain setup.

-   **Xw**, **Yw**, **Zw**: The wedge coordinate system.

    -   **origin**: Located in the middle and at the bottom of its front
        (blue axis system in the figure)**.**

-   **Xe**, **Ye**, **Ze**: The elements coordinate system.

    -   **origin**: Located in the middle of the first element surface
        (green axis system in the figure).

-   **X**, **Y**, **Z**: Global referential axis. It is independent of
    the acquisition and serves to position the data on the specimen.
    These are akin to real-world coordinates defined by the user for the
    specimen.

    -   **origin**: Arbitrary and stays the same across files for a
        given specimen. It is usually defined by the user on the
        specimen with some marking or physical reference in the specimen
        environment.

    -   **[NOTE]{.underline}**: Currently, with NDE format version
        3.1.0, the **X**, **Y**, and **Z** axis are not used nor defined
        in the NDE file, but to show inspection results in 3D, one would
        have to translate everything to this coordinate system.

-   **beams** axis: An axis used in the HDF5 dataset rather than a
    physical one where each element contains one beam\'s positions and
    parameters. It is used in PAUT scenarios when the beams do not fit
    well in a **U**, **V** grid. Giving the coordinate by beams
    simplifies their use in these scenarios.

-   **Ultrasound** axis**:** Time-based information sampled by an
    ultrasonic acquisition system and the specified probe/wedge
    configuration.

## Wedge:

-   The wedge **origin** is centered at the bottom of its front face
    (blue axis system in the figure).

-   The wedge **skew angle** is defined by the angle between the wedge
    and **U** axis.

-   Its starting position is given in relation to the **U** and **V**
    axis with the **uCoordinateOffset** and **vCoordinateOffset**
    located in the domain setup at **wedges\[0\].positioning**.

-   The width, height and length of the wedge are found at
    **wedges\[0\].angleBeamWedge**.

## Probes:

-   Position of the first element is given by; **\"primaryOffset\"**,
    **\"secondaryOffset\"** and **\"tertiaryOffset\"** found in
    **wedges\[0\].angleBeamWedge.mountingLocations**. These three
    offsets are given in relation to the wedge coordinate system where
    the primary, secondary, and tertiary offsets are on the **Yw**,
    **Xw**, and **Zw** axis respectively.

-   **Reference Position** (.Ref): Probes are referenced to other
    objects (wedge, specimen, etc.) through this point. The reference
    position is defined from the position of the first element when the
    probe skew = 0.

![Diagram, engineering drawing Description automatically
generated](C:/Dev/nde-format/documents/media/image2.png){width="8.172222222222222in"
height="5.155555555555556in"}

![Chart Description automatically
generated](C:/Dev/nde-format/documents/media/image3.png){width="3.580903324584427in"
height="3.4006944444444445in"}

Figure 2 - Raster example with the **U** and **V** axis.

![Chart Description automatically
generated](C:/Dev/nde-format/documents/media/image4.png){width="4.799212598425197in"
height="3.828905293088364in"}

![Chart Description automatically
generated](C:/Dev/nde-format/documents/media/image4.png){width="7.5in"
height="5.341666666666667in"}![A picture containing text, sky, several
Description automatically
generated](C:/Dev/nde-format/documents/media/image5.png){width="7.5in"
height="6.470833333333333in"}

![Diagram, engineering drawing Description automatically
generated](C:/Dev/nde-format/documents/media/image6.png){width="7.8284722222222225in"
height="4.298611111111111in"}
