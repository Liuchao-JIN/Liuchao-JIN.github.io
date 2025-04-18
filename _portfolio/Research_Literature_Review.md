---
title: "Research Literature Review"
excerpt: ""
collection: portfolio
---



Machine learning for hierarchical architecture design
======
* Li Beichen et al. 2024. Computational discovery of microstructured composites with optimal stiffness-toughness trade-offs. Science Advances.
  * 
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/li2024computational.pdf)\]
* Sun Xiaohao et al. 2024. Machine learning and sequential subdomain optimization for ultrafast inverse design of 4D-printed active composite structures. Journal of the Mechanics and Physics of Solids.
  * Method: machine learning (ML) and sequential subdomain optimization (SSO)
  * Main idea: Divide the part you want to design into several sections, then design the allocation of each section in sequence (you can go through all of them) to achieve the desired effect. <img src='/files/essay/sun2024machine_1.jpg'>
  * Advantage:
    * Fast <img src='/files/essay/sun2024machine_2.jpg'>
    * Can design varying numbers of lengthwise voxels <img src='/files/essay/sun2024machine_3.jpg'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/sun2024machine.pdf)\]




4D Printing
======
* Centrifugal multimaterial 3D printing of multifunctional heterogeneous objects
  * Keywords: Multimaterial DLP
  * Summary: There are growing demands for multimaterial three-dimensional (3D) printing to manufacture 3D object where voxels with different properties and functions are precisely arranged. Digital light processing (DLP) is a high-resolution fast-speed 3D printing technology suitable for various materials. However, multimaterial 3D printing is challenging for DLP as the current multimaterial switching methods require direct contact onto the printed part to remove residual resin. Here we report a DLP-based centrifugal multimaterial (CM) 3D printing method to generate large-volume heterogeneous 3D objects where composition, property and function are programmable at voxel scale. Centrifugal force enables non-contact, high-efficiency multimaterial switching, so that the CM 3D printer can print heterogenous 3D structures in large area (up to 180 mm × 130 mm) made of materials ranging from hydrogels to functional polymers, and even ceramics. Our CM 3D printing method exhibits excellent capability of fabricating digital materials, soft robots, and ceramic devices.
  * Note:
    * Mechanics, displacement field, tension-shear combination of digital printed tunable modulus part: The CM 3D printer also enables us to design and fabricate digital materials where the mechanical properties can be tuned by controlling the spatial distribution of the hard and soft voxels (Fig. 3d). By increasing the content of hard voxels from 0 to 100%, the modulus of the printed digital material raises from 0.8 MPa to 1 GPa (Fig. 3e, Supplementary Fig. 14). The capability of printing digital materials allows us to use only two base materials to design and fabricate one single part that exhibits multiple mechanical properties at different locations (Fig. 3f). We further apply this unique capability to a four-dimensional (4D) printing demonstration (Fig. 3g–j) where the palm and five fingers of a hand are formed with different digital materials, and a layer of hydrogel is printed on the top of the hand (Fig. 3g, h). After placing the hand into water for 1 h, the swelling of the hydrogel layer drives the five fingers to bend to different angles due to the different modulus (Fig. 3i). The hand finally makes a fist after being placed into water for 6 h (Fig. 3j).
      * We cannot control how much the hydrogel will deform after heating. Maybe not need for heating.
      * The size of the part can be 20 mm × 40 mm × 5 mm, the voxel size is 0.5 mm × 0.5 mm × 0.5 mm
  * <img src='/files/essay/cheng2022centrifugal.webp'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/cheng2022centrifugal.pdf)\] \[[Web](https://doi.org/10.1038/s41467-022-35622-6)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/cheng2022centrifugal.txt)\]
  * Author: Jianxiang Cheng, Rong Wang, Zechu Sun, Qingjiang Liu, Xiangnan He, Honggeng Li, Haitao Ye, Xingxin Yang, Xinfeng Wei, Zhenqing Li, Bingcong Jian, Weiwei Deng, Qi Ge
  * Year: 2022
  * Journal: Nature Communications
* Codesign of Biobased Cellulose-Filled Filaments and Mesostructures for 4D Printing Humidity Responsive Smart Structures
  * Keywords: fused filament fabrication, biobased polymers, hygromorphs, material programming, adaptive architecture
  * Summary: The article discusses a new approach to 4D printing of hygromorphic smart structures that can respond to relative humidity (RH). Hygromorphic structures can autonomously change their shape in response to environmental changes in RH, and their applications include adaptive shading elements and weather-responsive building envelopes. Four-dimensional (4D) printing, which involves 3D printing of structures that can change their shape over time in response to external stimuli, is a suitable method for developing such structures. However, current material limitations in terms of printability, responsiveness, and mechanical properties are major bottlenecks in achieving reliable and repeatable humidity-responsive actuation. The article proposes a codesign method for 4D printing hygromorphic structures through fused filament fabrication, which involves the development of cellulose-filled filaments with varying stiffness and hygroresponsiveness and designed mesoscale structuring in printed elements. The prototypes developed in the study can fully transform in conditions of 35–90% RH, which corresponds to naturally occurring shifts in RH in daily and seasonal weather cycles, and their motion is fast, fully reversible, and repeatable in numerous cycles. The study demonstrates the potential of using 4D printing and natural resources for the development of functional humidity-responsive smart structures.
  * <img src='/files/essay/tahouni2023codesign.jpg'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/tahouni2023codesign.pdf)\] \[[Web](https://www.liebertpub.com/doi/10.1089/3dp.2022.0061)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/tahouni2023codesign.txt)\]
  * Author: Yasaman Tahouni, Tiffany Cheng, Silvia Lajewski, Johannes Benz, Christian Bonten, Dylan Wood, Achim Menges
  * Year: 2023
  * Journal: 3D Printing and Additive Manufacturing
* Programming 3D curved mesosurfaces using microlattice designs
  * Keywords: Rational inverse design of 3D shapes
  * Summary: Cheng et al. developed an inverse design method to achieve complex three-dimensional (3D surfaces through a subset of 2D films that are bonded together. Analytic modeling and computations to inverse design the 2D patterns allow for control of the final porosity. A wide range of examples are provided, including changes in the sign of the curvature. These structures can be fabricated from silicon, metals, chitosan, and polymers.
  * <img src='/files/essay/cheng2023programming.jpg'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/cheng2023programming.pdf)\] \[[Web](https://www.science.org/doi/full/10.1126/science.adf3824)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/cheng2023programming.txt)\]
  * Author: Cheng Xu, Fan Zhichao, Yao Shenglian, Jin Tianqi, Lv Zengyao, Lan Yu, Bo Renheng, Chen Yitong, Zhang Fan, Shen Zhangming
  * Year: 2023
  * Journal: Science
* 4D Thermo-Responsive Smart hiPSC-CM Cardiac Construct for Myocardial Cell Therapy
  * Keywords: 4D printing, shape memory, nanostructure, myocardial regeneration, cellularized patch, minimally invasive
  * Summary: 4D fabrication techniques have been utilized for advanced biomedical therapeutics due to their ability to create dynamic constructs that can transform into desired shapes on demand. The internal structure of the human cardiovascular system is complex, where the contracting heart has a highly curved surface that changes shape with the heart’s dynamic beating motion. Hence, 4D architectures that adjust their shapes as required are a good candidate to readily deliver cardiac cells into the damaged heart and/or to serve as self-morphing tissue scaffolds/patches for healing cardiac diseases. In this proof-of-concept in vitro study, a two-in-one 4D smart cardiac construct that integrates the functions of minimally invasive cell vehicles and in situ tissue patches was developed for repairing damaged myocardial tissue.
  * <img src='/files/essay/hann20234d.jpg'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/hann20234d.pdf)\] \[[Web](https://www.tandfonline.com/doi/full/10.2147/IJN.S402855)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/hann20234d.txt)\]
  * Author: Sung Yun Hann, Haitao Cui, Timothy Esworthy, Lijie Grace Zhang
  * Year: 2023
  * Journal: International Journal of Nanomedicine
* Tough PEGgels by In Situ Phase Separation for 4D Printing
  * Keywords: 4D printing, gels materials
  * Summary: Polymer gels, consisting of cross-linked polymer network systems swollen by a solvent, show great potential in biomedicine, flexible electronics, and artificial muscles, due to their tissue-like mechanical properties. Due to the presence of a large amount of solvent, the improvement of the mechanical properties of the polymer gel is a challenge. Moreover, combining high toughness with useful properties, such as 3D printability or shape-memory, in one polymer gel system is even more challenging. In this study, a simple and efficient method is developed for the fabrication of tough polymer gels by polymerizing 2-hydroxyethyl methacrylate (HEMA) in a mixture of poly(ethylene glycol) (PEG) and poly(propylene glycol) (PPG). The polymerized elastic networkpresents distinct compatibility with PEG (compatible) and PPG (poorly compatible), resulting in in-situ phase separation at the microscale. The resulting phase-separated gel demonstrates high strength (8.0 MPa), favorable fracture strain (430%), and large toughness (17.0 MJ m−3). The separated hard phasewith a high glass transition temperature (75 °C) endows the whole soft polymer gel with the property of shape memory at room temperature. Finally, the fabrication of tunable tough PEGgels is combined with 3D printing as well as with shape memory properties, demonstrating the use of PEGgels for 4D printing.
  * <img src='/files/essay/wang2023tough.jpg'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/wang2023tough.pdf)\] \[[Web](https://onlinelibrary.wiley.com/doi/10.1002/adfm.202300947)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/wang2023tough.txt)\]
  * Author: Zhenwu Wang, Matthias Heck, Wenwu Yang, Manfred Wilhelm, Pavel A. Levkin
  * Year: 2023
  * Journal: Advanced Functional Materials
* 4D printing of multiple shape memory polymer and nanocomposites with biocompatible, programmable and selectively actuated properties
  * Keywords: 3 layer of SMP to control the shape change of materials continuously
  * Summary: 4D printing of poly (D,L-lactide-co-trimethylene carbonate) (PLMC)/poly (trimethylene carbonate) (PTMC)/Fe3O4 multi-material with multiple shape-changing capabilities under sequential stimuli of remotely magnetic field and heat was achieved. At first, we optimized the composition of pure SMP to fine tune the multiple shape memory effect and quantitatively characterized the shape recovery by stepwise heating. Then with the addition of Fe3O4 nanoparticles, the multi-material distribution of 4D printed structure consisting of multiple-SMP and its nanocomposites was designed. The integration of multi-material additive manufacturing with multiple shape memory effect extends the shape transformation to quintuple complex shapes with accurate and local controllability under selective multi-stimuli. The 4D printed multiple-SMP and its nanocomposites with simultaneously thermo- and magnetic- responsive shape-changing capability also demonstrated excellent biocompatibility. This work thus offers a feasible and robust approach for 4D printing of multi-functional devices for broad applications in entertainment, robotics, biomedical field and beyond.
  * <img src='/files/essay/wan20224d.jpg'> <img src='/files/essay/wan20224d_2.jpg'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/wan20224d.pdf)\] \[[Web](https://www.sciencedirect.com/science/article/pii/S221486042200094X?casa_token=BPlauJioENMAAAAA:h4llGSvKvKyO8rii8D6XjHXvunPXbNdZRnyLVK0GtIlc72yXaqw0hRF7nDj4EGGb5zvqR1Lf-do)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/wan20224d.txt)\]
  * Author: Xue Wan, Yang He, Yanju Liu, Jinsong Leng
  * Year: 2022
  * Journal: Additive Manufacturing
* 4D Printing of Hydrogels: A Review
  * Keywords: 4D Printing, Hydrogels, review
  * Summary: 3D printing permits the construction of objects by layer-by-layer deposition of material, resulting in precise control of the dimensions and properties of complex printed structures. Although 3D printing fabricates inanimate objects, the emerging technology of 4D printing allows for animated structures that change their shape, function, or properties over time when exposed to specific external stimuli after fabrication. Among the materials used in 4D printing, hydrogels have attracted growing interest due to the availability of various smart hydrogels. The reversible shape-morphing in 4D printed hydrogel structures is driven by a stress mismatch arising from the different swelling degrees in the parts of the structure upon application of a stimulus. This review provides the state-of-the-art of 4D printing of hydrogels from the materials perspective. First, the main 3D printing technologies employed are briefly depicted, and, for each one, the required physico-chemical properties of the precursor material. Then, the hydrogels that have been printed are described, including stimuli-responsive hydrogels, non-responsive hydrogels that are sensitive to solvent absorption/desorption, and multimaterial structures that are totally hydrogel-based. Finally, the current and future applications of this technology are presented, and the requisites and avenues of improvement in terms of material properties are discussed.
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/champeau20204d.pdf)\] \[[Web](https://onlinelibrary.wiley.com/doi/full/10.1002/adfm.201910606?casa_token=-IZ6jdDGeLoAAAAA%3Aoop7ndD0sXeg8nrU-9PpMrcMutFFkBPaiAReB-HVmhnw67MpyVX8PdtX9F7YYJGmLxl-pqd2tFHEkE_5Fw)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/champeau20204d.txt)\]
  * Author: Mathilde Champeau, Daniel Alves Heinze, Thiago Nunes Viana, Edcarlos Rodrigues de Souza, Anne Cristine Chinellato, Silvia Titotto
  * Year: 2020
  * Journal: Advanced Functional Materials
* Hydrophilic/Hydrophobic Composite Shape-Shifting Structures
  * Keywords: active structures solvent-responsive structures digital light processing 3D printing 4D printing
  * Summary: Swelling-induced shape transformation has been widely investigated and applied to the design and fabrication of smart polymer devices, such as soft robotics, biomedical devices, and origami patterns. Previous shape-shifting designs using soft hydrogels have several limitations, including relatively small actuation force, slow responsive speed, and relatively complicated fabrication process. In this paper, we develop a novel hydrophilic/hydrophobic composite structure by using photopolymers. The rubbery nature of the materials used in this composite provides desirable actuation speed and actuation force. The photocurable polymer system could be easily patterned by using the digital light processing technique. Experiments and theoretical analysis were conducted to study the actuation process. We also fabricated several three-dimensional water-responsive shape-shifting structures, including structures with sequential actuation behavior. Finally, the directional bending behavior of the hydrophilic/hydrophobic bilayer plate was investigated.
  * <img src='/files/essay/zhao2018hydrophilic.gif'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/zhao2018hydrophilic.pdf)\] \[[Web](https://pubs.acs.org/doi/10.1021/acsami.8b02444)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/zhao2018hydrophilic.txt)\]
  * Author: Zeang Zhao, Xiao Kuang, Chao Yuan, H. Jerry Qi, Daining Fang
  * Year: 2018
  * Journal: ACS applied materials & interfaces
* Responsive cellulose-hydrogel composite ink for 4D printing
  * Keywords: 4D materials, Additive manufacturing, Composite morphing, Stimuli-responsive, Cellulose-hydrogel
  * Summary: Sustainable and cost-effective solutions are crucial for the widespread adoption of 4D printing technology. This paper focuses on the development of a cellulose-hydrogel composite ink for additive manufacture, presenting the development and physical characterisation (stability, swelling potential and rheology) of the cellulose-hydrogel composite to establish its suitability for 4D printing of responsive structures. The use of a carboxymethyl cellulose (CMC) hydrocolloid with incorporated cellulose pulp fibres resulted in an ink with a high total cellulose content (fibre volume fraction ≈50% for the dehydrated composite) and good dispersion of fibres within the hydrogel matrix. The composite ink formulation developed in this study permitted smooth extrusion using an open source 3D printer to achieve controlled material placement in 3D space while retaining the functionality of the cellulose. The addition of montmorillonite clay not only resulted in enhanced storage stability of the composite ink formulations but also had a beneficial effect on the extrusion characteristics. The ability to precisely apply the ink via 3D printing was demonstrated through fabrication of a complex structure capable of morphing according to pre-determined design rules in response to hydration/dehydration.
  * <img src='/files/essay/mulakkal2018responsive.jpg'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/mulakkal2018responsive.pdf)\] \[[Web](https://www.sciencedirect.com/science/article/pii/S0264127518307032?pes=vor)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/mulakkal2018responsive.txt)\]
  * Author: Manu C. Mulakkal, Richard S. Trask, Valeska P. Ting, Annela M. Seddon
  * Year: 2018
  * Journal: Materials & Design
* Biomimetic 4D printing
  * Keywords: hydrogel 4D printing
  * Summary: Shape-morphing systems can be found in many areas, including smart textiles, autonomous robotics, biomedical devices, drug delivery and tissue engineering. The natural analogues of such systems are exemplified by nastic plant motions, where a variety of organs such as tendrils, bracts, leaves and flowers respond to environmental stimuli (such as humidity, light or touch) by varying internal turgor, which leads to dynamic conformations governed by the tissue composition and microstructural anisotropy of cell walls. Inspired by these botanical systems, we printed composite hydrogel architectures that are encoded with localized, anisotropic swelling behaviour controlled by the alignment of cellulose fibrils along prescribed four-dimensional printing pathways. When combined with a minimal theoretical framework that allows us to solve the inverse problem of designing the alignment patterns for prescribed target shapes, we can programmably fabricate plant-inspired architectures that change shape on immersion in water, yielding complex three-dimensional morphologies.
  * <img src='/files/essay/sydney2016biomimetic.webp'> <img src='/files/essay/sydney2016biomimetic_1.webp'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/sydney2016biomimetic.pdf)\] \[[Web](https://www.nature.com/articles/nmat4544)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/sydney2016biomimetic.txt)\]
  * Author: A. Sydney Gladman, Elisabetta A. Matsumoto, Ralph G. Nuzzo, L. Mahadevan, Jennifer A. Lewis
  * Year: 2016
  * Journal: Nature materials
* Active Printed Materials for Complex Self-Evolving Deformations
  * Keywords: hydrogel 4D printing
  * Summary: We propose a new design of complex self-evolving structures that vary over time due to environmental interaction. In conventional 3D printing systems, materials are meant to be stable rather than active and fabricated models are designed and printed as static objects. Here, we introduce a novel approach for simulating and fabricating self-evolving structures that transform into a predetermined shape, changing property and function after fabrication. The new locally coordinated bending primitives combine into a single system, allowing for a global deformation which can stretch, fold and bend given environmental stimulus.
  * <img src='/files/essay/raviv2014active.webp'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/raviv2014active.pdf)\] \[[Web](https://www.nature.com/articles/srep07422)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/raviv2014active.txt)\]
  * Author: Dan Raviv, Wei Zhao, Carrie McKnelly, Athina Papadopoulou, Achuta Kadambi, Boxin Shi, Shai Hirsch, Daniel Dikovsky, Michael Zyracki, Carlos Olguin, Ramesh Raskar, Skylar Tibbits
  * Year: 2014
  * Journal: Scientific reports
* Multimaterial 4D Printing with Tailorable Shape Memory Polymers
  * Keywords: Multimaterial 4D Printing
  * Summary: We present a new 4D printing approach that can create high resolution (up to a few microns), multimaterial shape memory polymer (SMP) architectures. The approach is based on high resolution projection microstereolithography (PμSL) and uses a family of photo-curable methacrylate based copolymer networks. We designed the constituents and compositions to exhibit desired thermomechanical behavior (including rubbery modulus, glass transition temperature and failure strain which is more than 300% and larger than any existing printable materials) to enable controlled shape memory behavior. We used a high resolution, high contrast digital micro display to ensure high resolution of photo-curing methacrylate based SMPs that requires higher exposure energy than more common acrylate based polymers. An automated material exchange process enables the manufacture of 3D composite architectures from multiple photo-curable SMPs. In order to understand the behavior of the 3D composite microarchitectures, we carry out high fidelity computational simulations of their complex nonlinear, time-dependent behavior and study important design considerations including local deformation, shape fixity and free recovery rate. Simulations are in good agreement with experiments for a series of single and multimaterial components and can be used to facilitate the design of SMP 3D structures.
  * <img src='/files/essay/ge2016multimaterial.webp'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/ge2016multimaterial.pdf)\] \[[Web](https://www.nature.com/articles/srep31110)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/ge2016multimaterial.txt)\]
  * Author: Qi Ge, Amir Hosein Sakhaei, Howon Lee, Conner K. Dunn, Nicholas X. Fang, Martin L. Dunn
  * Year: 2016
  * Journal: Scientific Reports
* Advances in 4D printed shape memory composites and structures: Actuation and application
  * Keywords: 4D printing, shape memory composites
  * Summary: Shape memory polymer composites (SMPCs) are a type of smart material that can change shapes under the stimulation of the external environment, and they have great potential in aerospace, biomedical, robotics, and electronic devices due to their advantages of high strength and toughness, lightweight, impact resistance, corrosion resistance, and aging resistance. 4D printing technology has provided new opportunities for the further development of smart materials. The addition of various fillers enriches the variety of printable materials and provides composites with different properties and functions. The combination of SMPCs and printing technologies realizes the structure-function integration. This paper introduces the emergence and development of 4D printing technologies, the preparation methods and properties of SMPCs for 4D printing; as well as the research progress and potential application of 4D printable SMPCs in recent years in terms of thermal, electrical, magnetic, and optical driving. Finally, the existing problems and future development of 4D printable SMPCs are discussed.
  * <img src='/files/essay/wang2023advances.png'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/wang2023advances.pdf)\] \[[Web](https://link.springer.com/article/10.1007/s11431-022-2255-0)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/wang2023advances.txt)\]
  * Author: LinLin Wang, FengHua Zhang, ShanYi Du, JinSong Leng
  * Year: 2023
  * Journal: Science China Technological Sciences
* 4D Printing: A Review on Recent Progresses
  * Keywords: four-dimensional (4D) printing, additive manufacturing, smart materials, shape memory polymer
  * Summary: Since the late 1980s, additive manufacturing (AM), commonly known as three-dimensional (3D) printing, has been gradually popularized. However, the microstructures fabricated using 3D printing is static. To overcome this challenge, four-dimensional (4D) printing which defined as fabricating a complex spontaneous structure that changes with time respond in an intended manner to external stimuli. 4D printing originates in 3D printing, but beyond 3D printing. Although 4D printing is mainly based on 3D printing and become an branch of additive manufacturing, the fabricated objects are no longer static and can be transformed into complex structures by changing the size, shape, property and functionality under external stimuli, which makes 3D printing alive. Herein, recent major progresses in 4D printing are reviewed, including AM technologies for 4D printing, stimulation method, materials and applications. In addition, the current challenges and future prospects of 4D printing were highlighted.
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/chu20204d.pdf)\] \[[Web](https://www.mdpi.com/2072-666X/11/9/796)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/chu20204d.txt)\]
  * Author: Honghui Chu, Wenguang Yang, Lujing Sun, Shuxiang Cai, Rendi Yang, Wenfeng Liang, Haibo Yu, Lianqing Liu
  * Year: 2020
  * Journal: Micromachines
* 4D Printing in Pharmaceutics and Biomedical Applications
  * Keywords: 4D Printing, Pharmaceutics, Biomedical
  * Summary: 3D printing (3DP) has made significant advancements in the past decade in the fabrication of complex objects that are based on biomaterials. Although 3D-printed constructs were promising for biomedical applications, they fell short due to their inability to accurately mimic dynamic human tissues. 4D printing (4DP) is a breakthrough delivery system that integrates “time” into the conventional concept of 3DP to address the dynamic healing and regeneration of human tissues. In that way, additive manufacturing (AM) goes from 3DP to 4DP and implicates the use of stimuli-responsive materials. With its ability to create a wide range of useful biomedical products, 4DP has become an important tool in biomedical engineering. The purpose of this chapter is to present the concept of 4D bioprinting and the recent developments in smart materials, which can be actuated by different stimuli and can be used to develop biomimicry materials and structures with significant implications for pharmaceutics and biomedical research, as well as perspectives for the future.
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/naniz20234d.pdf)\] \[[Web](https://link.springer.com/chapter/10.1007/978-3-031-26908-0_9)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/naniz20234d.txt)\]
  * Author: Moqaddaseh Afzali Naniz, Mohsen Askari, Ali Zolfagharian, Mahdi Bodaghi
  * Year: 2023
  * Book: Nano-and Microfabrication Techniques in Drug Delivery: Recent Developments and Future Prospects
* Flexible, biocompatible and highly conductive MXene-graphene oxide film for smart actuator and humidity sensor
  * Keywords: MXene, Graphene oxide, Actuator, Moisture gradient, Monitoring respiration, Humidity sensor
  * Summary: The evaporation of water occurs ubiquitously on earth. Hence, smart materials that can directly convert signals generated via water stimulation into mechanical motion have attracted wide attention. However, it is still a challenge to develop novel functional materials with fast response, large scale deformation, and long-term stability for moisture-gradient actuators. Here, a flexible, conductive, layer-structured homogenous Ti3C2TX MXene-graphene oxide (MGO) film-based moisture-driven actuator and humidity sensor were fabricated. The oxygen groups and d-spacing could be effectively adjusted by MXene/GO composition ratio, thereby tuning the actuation performance. MGO3 (MXene/GO = 3) displayed a large bending angle, and reversible deformation. And the bending speed of MGO3 is up to 32°s−1. Furthermore, MGO3 actuation displayed long-term stability via suppression of MXene oxidation by the introduction of GO and showed good cycling stability. MGO3 actuators are constructed, which could mimic the blooming of flower, lifting and carrying objects, and be used as a non-contact control switch. In addition, MGO3 showed a linear sensitive response to humidity and excellent biocompatibility which make it suitable for respiratory monitoring. This work demonstrated that flexible, biocompatibility and conductive MGO films have broad application prospects in the fields of smart actuators, sensing devices, and biology and health care.
  * <img src='/files/essay/jia2021flexible.jpg'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/jia2021flexible.pdf)\] \[[Web]()\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/jia2021flexible.txt)\]
  * Author: Guangwen Jia, Ao Zheng, Xiao Wang, Lu Zhang, Ling Li, Chenxing Li, Yan Zhang, Lingyan Cao
  * Year: 2021
  * Journal: Sensors and Actuators B: Chemical
* 4D printing soft robotics for biomedical applications
  * Keywords: 4D printing, soft robotics, biomedical applications, review
  * Summary: Soft robotics has grown rapidly as an attractive manufacturing technique at micro/nanoscales, especially in the field of biomedical engineering, due to its mobility and compact size. As an emerging additive manufacturing technique, four-dimensional (4D) printing can replicate natural physio-mechanical changes over time leading the transition from static to dynamic. As such, 4D printing is widely investigated and applied in fields ranging from mechanical engineering and material science, to biomedical engineering. By combining the unique ability of 4D printing to create dynamic morphological changes under certain stimuli and biocompatible soft material based micro-/nanorobots, a new promising platform with precise controllability and unlimited reversible actuation is expected to enhance the role of 4D printing for biomedical applications. In this review, we systematically introduce soft robots and further summarize current 4D printing approaches to fabricate soft robots. Moreover, a broad scope of potential applications of 4D soft robots in biomedical engineering is exclusively discussed. We finally conclude with current challenges and limitations, as well as future directions, for 4D printing soft robot technology.
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/hann20204d.pdf)\] \[[Web](https://doi.org/10.1016/j.addma.2020.101567)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/hann20204d.txt)\]
  * Author: Sung Yun Hann, Haitao Cui, Margaret Nowicki, Lijie Grace Zhang
  * Year: 2020
  * Journal: Additive Manufacturing
* Programming a crystalline shape memory polymer network with thermo- and photo-reversible bonds toward a single-component soft robot
  * Keywords: 4D printing, soft robotics
  * Summary: The need to support the two most basic functions [three-dimensional (3D)–shaped support and actuation] independently for a typical robot demands that at least two components should be used in its construction. Therefore, component assembly is unavoidable despite the ultimate dream of creating assembly-free robots. We devise a strategy that uses a programmable crystalline shape memory polymer with thermo- and photo-reversible bonds to create a single-component robot. The global 3D-shaped structural support is fabricated via a plasticity-based origami technique enabled by the thermo-reversible bonds. More critically, **precisely controlled localized actuation can be programmed into the 3D origami via spatially defined reversible shape memory using the photo-reversible bonds. The overall result is that a polymer thin film can be programmed into various soft robots including a 3D crane and an elephant**. Besides reversible shape memory, other types of actuation mechanisms can be potentially introduced via a similar principle. Thus, our strategy represents a general method to create single-component soft robots.
  * <img src='/files/essay/jin2018programming.jpg'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/jin2018programming.pdf)\] \[[Web](https://doi.org/10.1126/sciadv.aao3865)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/jin2018programming.txt)\]
  * Author: Binjie Jin, Huijie Song, Ruiqi Jiang, Jizhou Song, Qian Zhao, Tao Xie
  * Year: 2018
  * Journal: Science Advances
* Closed-loop 4D-printed soft robots
  * Keywords: 4D printing, soft robotics, review
  * Summary: Soft robotics is a recent and rapidly growing field of research that encompasses emerging advances in functional materials, fabrication, modelling, and performance control with important applications in the manipulation of fragile objects. The incorporation of three-dimensional (3D) printing into the fabrication of soft robots facilitates the customisation of their functions by the strategic placement of functional materials into locations that may be inaccessible by conventional manufacturing methods. Most current 3D-printed soft robots, however, are fabricated without considering their autonomy. In four-dimensional (4D) printing, the functionality of a soft robot is introduced during the printing process. To control this functionality, specific functional materials are embedded in desired locations. Four-dimensional printing and machine learning techniques provide new possibilities for developing stand-alone closed-loop 4D-printed soft robots. This review paper presents the current approaches employed to design and construct 4D-printed soft robots with customised geometrical, functional, and control properties.
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/zolfagharian2020closed.pdf)\] \[[Web](https://doi.org/10.1016/j.matdes.2019.108411)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/Essayzolfagharian2020closed/.txt)\]
  * Author: Ali Zolfagharian, Akif Kaynak, Abbas Kouzani
  * Year: 2020
  * Journal: Materials & Design
* Shapeshifting: Reversible Shape Memory in Semicrystalline Elastomers
  * Keywords: **Reversible Shape Memory**
  * Summary: We present a general strategy for enabling reversible shape transformation in semicrystalline shape memory (SM) materials, which integrates three different SM behaviors: conventional one-way SM, two-way reversible SM, and one-way reversible SM. While two-way reversible shape memory (RSM) is observed upon heating and cooling cycles, the one-way RSM occurs upon heating only. Shape reversibility is achieved through partial melting of a crystalline scaffold which secures memory of a temporary shape by leaving a latent template for recrystallization. This behavior is neither mechanically nor structurally constrained, thereby allowing for multiple switching between encoded shapes without applying any external force, which was demonstrated for different shapes including hairpin, coil, origami, and a robotic gripper. Fraction of reversible strain increases with cross-linking density, reaching a maximum of ca. 70%, and then decreases at higher cross-linking densities. This behavior has been shown to correlate with efficiency of securing the temporary shape.
  * <img src='/files/essay/zhou2014shapeshifting.png'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/zhou2014shapeshifting.pdf)\] \[[Web](https://doi.org/10.1021/ma4023185)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/zhou2014shapeshifting.txt)\]
  * Author: Jing Zhou, Sara A. Turner, Sarah M. Brosnan, Qiaoxi Li, Jan-Michael Y. Carrillo, Dmytro Nykypanchuk, Oleg Gang, Valerie S. Ashby, Andrey V. Dobrynin, Sergei S. Sheiko
  * Year: 2014
  * Journal: Macromolecules
* Switchable Micropatterned Surface Topographies Mediated by Reversible Shape Memory
  * Keywords: **Reversible Shape Memory**
  * Summary: Reversibly switching topography on micrometer length scales greatly expands the functionality of stimuli-responsive substrates. Here we report the first usage of reversible shape memory for the actuation of two-way transitions between microscopically patterned substrates, resulting in corresponding modulations of the wetting properties. Reversible switching of the surface topography is achieved through partial melting and recrystallization of a semi-crystalline polyester embossed with microscopic features. This behavior is monitored with atomic force microscopy (AFM) and contact angle measurements. We demonstrate that the magnitude of the contact angle variations depends on the embossment pattern.
  * <img src='/files/essay/turner2014switchable.jpg'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/turner2014switchable.pdf)\] \[[Web](https://doi.org/10.1021/am501970d)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/turner2014switchable.txt)\]
  * Author: Sara A. Turner, Jing Zhou, Sergei S. Sheiko, Valerie Sheares Ashby
  * Year: 2014
  * Journal: Macromolecules
* Photothermally and magnetically controlled reconfiguration of polymer composites for soft robotics
  * Keywords: 4D printing, soft robotics, reversible
  * Summary: New materials are advancing the field of soft robotics. Composite films of magnetic iron microparticles dispersed in a shape memory polymer matrix are demonstrated for reconfigurable, remotely actuated soft robots. The composite films simultaneously respond to magnetic fields and light. Temporary shapes obtained through combined magnetic actuation and photothermal heating can be locked by switching off the light and magnetic field. Subsequent illumination in the absence of the magnetic field drives recovery of the permanent shape. In cantilevers and flowers, multiple cycles of locking and unlocking are demonstrated. Scrolls show that the permanent shape of the film can be programmed, and they can be frozen in intermediate configurations. Bistable snappers can be magnetically and optically actuated, as well as biased, by controlling the permanent shape. Grabbers can pick up and release objects repeatedly. Simulations of combined photothermal heating and magnetic actuation are useful for guiding the design of new devices.
  * <img src='/files/essay/miao20164d.jpg'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/miao20164d.pdf)\] \[[Web](https://doi.org/10.1126/sciadv.aaw2897)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/miao20164d.txt)\]
  * Author: Jessica A.-C. Liu, Jonathan H. Gillen, Sumeet R. Mishra, Benjamin A. Evans, Joseph B. Tracy
  * Year: 2019
  * Journal: Science Advances
* Computational design for 4D printing of topology optimized multi-material active composites
  * Keywords: topology optimization for 4D printing
  * Summary: Recent efforts on design for four-dimensional (4D) printing have considered the spatial arrangement of smart materials and energy stimuli. The development of multifunctional structures and their desired mechanical/actuation performances require tackling 4D printing from a multi-material design perspective. With the materials distributions there is an opportunity to increase the spectrum of design concepts with computational approaches. The main goal being to achieve the “best” distribution of material properties in a voxelized structure, a computational framework that consists of a finite element analysis-based evolutionary algorithm is presented. It fuses the advantages of optimizing both the materials distribution and material layout within a design space via topology optimization to solve the inverse design problem of finding an optimal design to achieve a target shape change by integrating void voxels. The results demonstrate the efficacy of the proposed method in providing a highly capable tool for the design of 4D-printed active composites.
  * <img src='/files/essay/athinarayanarao2023computational.webp'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/athinarayanarao2023computational.pdf)\] \[[Web](https://doi.org/10.1038/s41524-022-00962-w)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/athinarayanarao2023computational.txt)\]
  * Author: Darshan Athinarayanarao, Romaric Prod’hon, Dominique Chamoret, H. Jerry Qi, Mahdi Bodaghi, Jean-Claude André, Frédéric Demoly
  * Year: 2023
  * Journal: npj Computational Materials

<!-- *
  * Keywords:
  * Summary:
  * <img src='/files/essay/.'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/.pdf)\] \[[Web](https://doi.org/)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/.txt)\]
  * Author:
  * Year:
  * Journal:  -->

3D Printing
======
* Rotational multimaterial printing of filaments with subvoxel control
  * Keywords: origami robots by embedding sensing, computing, and actuating in compliant, conductive materials, without requiring semiconductor-based electronics
  * Summary: Larson et al. developed a rotational multimaterial 3D printing method for creating helical filaments. Using this new approach, the team designed and fabricated artificial muscles and springy lattices for use in soft robotics and structural applications.
  * <img src='/files/essay/larson2023rotational.webp'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/larson2023rotational.pdf)\] \[[Web](https://www.nature.com/articles/s41586-022-05490-7?utm_medium=organic_social&utm_source=wechat&utm_campaign=CONR_PF020_BAWG_AP_CNCM_002E8_natvideo)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/larson2023rotational.txt)\]
  * Author: Natalie M. Larson, Jochen Mueller, Alex Chortos, Zoey S. Davidson, David R. Clarke, Jennifer A. Lewis
  * Year: 2023
  * Journal: Nature Communications
* Green 3D-printed lattice-shaped suspension arms for RC cars
  * Keywords: Suspension system, Sustainable design, Lattice structure, Modal analysis, 3D printing, Fused deposition modelling
  * Summary: This study proposes sustainable design recommendations for lattice-shaped plastic suspension arms in remote-controlled automobiles to achieve high performance while minimizing waste. Solid, re-entrant honeycomb, face-centred cubic, hexagonal honeycomb, hexagonal prism diamond, simple cubic, triangular honeycomb, and lattice from a volume mesh are implemented on suspension arms to investigate their effects on flexural strength and frequency response. Finite element analysis and experiments are conducted to examine the deformation and vibration frequency. The results show the relationships between mass, stiffness, and vibration frequency.
  * <img src='/files/essay/lalegani2023green.webp'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/lalegani2023green.pdf)\] \[[Web](https://scholar.googleusercontent.com/scholar.bib?q=info:1Hs4alhNO8sJ:scholar.google.com/&output=citation&scisdr=Cpv3NHohEOusl65l2yU:AJ9-iYsAAAAAZDdjwyXDg5uWy5zmkY1CxQRA64k&scisig=AJ9-iYsAAAAAZDdjw_CUfZBBfeVxKzymSaTZes0&scisf=4&ct=citation&cd=0&hl=en)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/lalegani2023green.txt)\]
  * Author: Mohammadreza Lalegani Dezaki, Mahdi Bodaghi, Ahmad Serjouei, Ali Zolfagharian
  * Year: 2023
  * Journal: Progress in Additive Manufacturing
* Tough, Transparent, 3D-Printable, and Self-Healing Poly(ethylene glycol)-Gel (PEGgel)
  * Keywords: 3D printing, soft robots
  * Summary: Polymer gels, such as hydrogels, have been widely used in biomedical applications, flexible electronics, and soft machines. Polymer network design and its contribution to the performance of gels has been extensively studied. In this study, the critical influence of the solvent nature on the mechanical properties and performance of soft polymer gels is demonstrated. A polymer gel platform based on poly(ethylene glycol) (PEG) as solvent is reported (PEGgel). Compared to the corresponding hydrogel or ethylene glycol gel, the PEGgel with physically cross-linked poly(hydroxyethyl methacrylate-co-acrylic acid) demonstrates high stretchability and toughness, rapid self-healing, and long-term stability. Depending on the molecular weight and fraction of PEG, the tensile strength of the PEGgels varies from 0.22 to 41.3 MPa, fracture strain from 12% to 4336%, modulus from 0.08 to 352 MPa, and toughness from 2.89 to 56.23 MJ m–3. Finally, rapid self-healing of the PEGgel is demonstrated and a self-healing pneumatic actuator is fabricated by 3D-printing. The enhanced mechanical properties of the PEGgel system may be extended to other polymer networks (both chemically and physically cross-linked). Such a simple 3D-printable, self-healing, and tough soft material holds promise for broad applications in wearable electronics, soft actuators and robotics.
  * <img src='/files/essay/wang2022tough_1.webp'> <img src='/files/essay/wang2022tough_2.jpg'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/wang2022tough.pdf)\] \[[Web](https://onlinelibrary.wiley.com/doi/full/10.1002/adma.202107791)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/wang2022tough.txt)\]
  * Author: Zhenwu Wang, Haijun Cui, Modan Liu, Stephan L. Grage, Maxi Hoffmann, Elaheh Sedghamiz, Wolfgang Wenzel, Pavel A. Levkin
  * Year: 2022
  * Journal: Advanced Materials
* Rheological behaviour of different composite materials for additive manufacturing of 3D bone scaffolds
  * Keywords: Additive manufacturing, Polymer-ceramic blends, Printability, Rheology
  * Summary: The aim of the paper is the investigation of rheological behaviour of polymer and composite blends regularly used for the production of scaffolds for bone tissue applications with the use of additive manufacturing. Poly-ε-caprolactone (PCL), hydroxyapatite (HA), β-tri-calcium phosphate (TCP) and Bioglass 45S5 blends containing different ceramic concentrations (10wt%, 15wt% and 20wt%) were prepared with the use of melt blending procedure and investigated with the use of oscillation and rotational rheology tests. Results are showing that all blends are presenting viscoelastic behaviour with higher viscous modulus, compared with elastic modulus for low frequencies, with this difference reducing while the frequency is increasing. All blends are presenting shear-thinning behaviour suitable for use with additive manufacturing methods. Viscous and elastic modulus are increasing by adding ceramic particles. Results are presenting that PCL/HA blends of the same material concentration are presenting higher elastic modulus properties compared with the other blends, while PCL/Bioglass blends are presenting lower loss factor, lower relaxation time and lower shear viscosity making them easier to handle during the printing procedure.
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/daskalakis2023rheological.pdf)\] \[[Web](https://www.sciencedirect.com/science/article/pii/S2238785423006944)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/daskalakis2023rheological.txt)\]
  * Author: Evangelos Daskalakis, Mohamed H. Hassan, Abdalla M. Omar, Glen Cooper, Andrew Weightman, Paulo Bartolo
  * Year: 2023
  * Journal: Journal of Materials Research and Technology
* MechSense: A Design and Fabrication Pipeline for Integrating Rotary Encoders into 3D Printed Mechanisms
  * Keywords: 3D printed mechanisms, printed electronics, capacitive sensing.
  * Summary: Researchers, including Associate Professor Stephanie Mueller, have created a system that enables makers to incorporate sensors directly into rotational mechanisms with only one pass in a 3D printer. This system gives rotational tools like the gears inside a motor the ability to sense their angular position, rotation speed, and direction of rotation.
  * <img src='/files/essay/mueller2023mechsense.jpg'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/alalawi2023mechsense.pdf)\] \[[Web](https://dl.acm.org/doi/abs/10.1145/3544548.3581361)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/alalawi2023mechsense.txt)\]
  * Author: Marwa Alalawi, Noah Pacik-Nelson, Junyi Zhu, Ben Greenspan, Andrew Doan, Brandon M Wong, Benjamin Owen-Block, Shanti Kaylene Mickens, Wilhelm Jacobus Schoeman, Michael Wessely, Andreea Danielescu, Stefanie Mueller
  * Year: 2023
* Multimaterial magnetically assisted 3D printing of composite materials
  * Keywords: DIW, make printing part have magnetical direction
  * Summary: 3D printing has become commonplace for the manufacturing of objects with unusual geometries. Recent developments that enabled printing of multiple materials indicate that the technology can potentially offer a much wider design space beyond unusual shaping. Here we show that a new dimension in this design space can be exploited through the control of the orientation of anisotropic particles used as building blocks during a direct ink-writing process. Particle orientation control is demonstrated by applying low magnetic fields on deposited inks pre-loaded with magnetized stiff platelets. Multimaterial dispensers and a two-component mixing unit provide additional control over the local composition of the printed material. The five-dimensional design space covered by the proposed multimaterial magnetically assisted 3D printing platform (MM-3D printing) opens the way towards the manufacturing of functional heterogeneous materials with exquisite microstructural features thus far only accessible by biological materials grown in nature.
  * <img src='/files/essay/kokkinis2015multimaterial.webp'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/kokkinis2015multimaterial.pdf)\] \[[Web](https://www.nature.com/articles/ncomms9643)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/kokkinis2015multimaterial.txt)\]
  * Author: Dimitri Kokkinis, Manuel Schaffner, André R. Studart
  * Year: 2015
  * Journal: Nature Communications
* A self-healing nanocomposite double network bacterial nanocellulose/gelatin hydrogel for three dimensional printing
  * Keywords: self-healing nanocomposite, 3D printing hydrogel
  * Summary: Extrusion-based three-dimensional (3D) printing of gelatin is important for additive manufactured tissue engineering scaffolds, but gelatin's thermal instability has remained an ongoing challenge. The gelatin tends to suddenly collapse at mild temperatures, which is a significant limitation for using it at physiological temperature of 37 °C. Hence, fabrication of a thermo-processable gelatin hydrogel adapted for extrusion-based additive manufacturing is still a challenge. To achieve this, a self-healing nanocomposite double-network (ncDN) gelatin hydrogel was fabricated with high thermo-processability, shear-thinning, mechanical strength, self-healing, self-recovery, and biocompatibility. To do this, amino group-rich gelatin was first created by combining gelatin with carboxyl methyl chitosan. Afterwards, a self-healing ncDN gelatin hydrogel was formed via an in-situ formation of imine bonds between the blend of gelatin/carboxyl methyl chitosan (Gel/CMCh) and dialdehyde-functionalized bacterial nanocellulose (dBNC). dBNC plays as nanofiber cross-linkers capable of simultaneously crosslinking and reinforcing the double networks of Gel/CMCh through formation of dynamic 3D imine bonds. Based on our findings, our self-healing ncDA gelatin hydrogel displayed great potential as a promising ink for additive manufactured tissue engineering scaffolds.
  * <img src='/files/essay/heidarian2023self.jpg'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/heidarian2023self.pdf)\] \[[Web](https://www.sciencedirect.com/science/article/pii/S0144861723003442)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/heidarian2023self.txt)\]
  * Author: Pejman Heidarian, Abbas Z. Kouzani
  * Year: 2023
  * Journal: Carbohydrate Polymers
* Scalable submicrometer additive manufacturing
  * Keywords: two-photon lithography
  * Summary: Saha et al. optimize a new parallel printing methodology that relies on ultrafast lasers. They show the ability to dramatically increase the speed of printing while maintaining submicrometer resolution.
  * <img src='/files/essay/saha2019scalable.jpg'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/saha2019scalable.pdf)\] \[[Web](https://www.science.org/doi/full/10.1126/science.aax8760)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/saha2019scalable.txt)\]
  * Author: Sourabh K. Saha, Dien Wang, Vu H. Nguyen, Yina Chang, James S. Oakdale, Shih-Chi Chen
  * Year: 2019
  * Journal: Science
* Three-dimensional nanofabrication via ultrafast laser patterning and kinetically regulated material assembly
  * Keywords: three-dimensional nanofabrication
  * Summary: Han et al. synthesized very finely detailed objects from a wide range of materials using femtosecond light sheets and nanoparticle-laden hydrogels. The strategy works for ceramics, polymers, metals, semiconductors, and other materials while still maintaining fine feature sizes. This technique could enable nanofabrication across different classes of materials.
  * <img src='/files/essay/han2022three.jpg'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/han2022three.pdf)\] \[[Web](https://www.science.org/doi/full/10.1126/science.abm8420)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/han2022three.txt)\]
  * Author: Fei Han, Songyun Gu1, Aleks Klimas, Ni Zhao, Yongxin Zhao, Shih-Chi Chen
  * Year: 2022
  * Journal: Science
* Ultrafast 3D nanofabrication via digital holography
  * Keywords: three-dimensional nanofabrication
  * Summary: There has been a compelling demand of fabricating high-resolution complex three-dimensional (3D) structures in nanotechnology. While two-photon lithography (TPL) largely satisfies the need since its introduction, its low writing speed and high cost make it impractical for many large-scale applications. We report a digital holography-based TPL platform that realizes parallel printing with up to 2000 individually programmable laser foci to fabricate complex 3D structures with 90 nm resolution. This effectively improves the fabrication rate to 2,000,000 voxels/sec. The promising result is enabled by the polymerization kinetics under a low-repetition-rate regenerative laser amplifier, where the smallest features are defined via a single laser pulse at 1 kHz. We have fabricated large-scale metastructures and optical devices of up to centimeter-scale to validate the predicted writing speed, resolution, and cost. The results confirm our method provides an effective solution for scaling up TPL for applications beyond laboratory prototyping.
  * <img src='/files/essay/ouyang2023ultrafast.webp'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/ouyang2023ultrafast.pdf)\] \[[Web](https://www.nature.com/articles/s41467-023-37163-y)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/ouyang2023ultrafast.txt)\]
  * Author: Wenqi Ouyang, Xiayi Xu, Wanping Lu, Ni Zhao, Fei Han, Shih-Chi Chen
  * Year: 2023
  * Journal: Nature Communications
* Ultrastrong and damage-tolerant ceramic architectures via 3D printing
  * Keywords: Damage tolerance, Ceramic architectures, 3D printing, Mechanical property, In situ compression
  * Summary: Ceramic materials have high mechanical strength and exceptional environmental stability, but are suboptimal for structural applications due to their inherent brittleness and low damage tolerance. Here, we report ultrastrong and damage-tolerant ceramic architectures that are designed based on Schwarz Primitive structures and manufactured by digital light processing (DLP)-based 3D printing. Through micro-computed tomography imaging and in situ compression experiments, we reveal that the step effect of 3D printing plays a role in crack initiation and propagation of the 3D-printed architectures. When the loading direction is perpendicular to the printing direction, the steps can induce cracks to propagate along the loading direction, which is beneficial to localize the cracks in the outer part of structure and enhance the structural strength and damage tolerance. After optimizing structural design and heat treatment process, the printed ceramic architecture can achieve compressive strength as high as 710 MPa at a relative density of 57.58 %. More importantly, the printed ceramic architecture exhibits excellent damage tolerance. It can bear more 20 cycles of 6 % compressive strain when 28 % of the structures have been damaged. The ceramic architecture can still bear the load without failure even when the degree of damage reaches 44 %. The superior mechanical properties make them have great potential in engineering applications that require high mechanical reliability.
  * <img src='/files/essay/wang2023ultrastrong.jpg'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/wang2023ultrastrong.pdf)\] \[[Web](https://www.sciencedirect.com/science/article/pii/S2214860422007503)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/wang2023ultrastrong.txt)\]
  * Author: Rong Wang, Haitao Ye, Jianxiang Cheng, Honggeng Li, Pengfei Zhu, Bo Li, Rong Fan, Juzheng Chen, Yang Lu, Qi Ge
  * Year: 2023
  * Journal: Additive Manufacturing
* Multimaterial Three-Dimensional Printing of Ultraviolet-Curable Ionic Conductive Elastomers with Diverse Polymers for Multifunctional Flexible Electronics
  * Keywords: flexible electronics ionic conductive elastomers multimaterial 3D printing digital light processing 4D printing
  * Summary: Ionic conductive elastomers (ICEs) are emerging stretchable and ionic conductive materials that are solvent-free and thus demonstrate excellent thermal stability. Three-dimensional (3D) printing that creates complex 3D structures in free forms is considered as an ideal approach to manufacture sophisticated ICE-based devices. However, the current technologies constrain 3D printed ICE structures in a single material, which greatly limits functionality and performance of ICE-based devices and machines. Here, we report a digital light processing (DLP)-based multimaterial 3D printing capability to seemly integrate ultraviolet-curable ICE (UV-ICE) with nonconductive materials to create ionic flexible electronic devices in 3D forms with enhanced performance. This unique capability allows us to readily manufacture various 3D flexible electronic devices. To demonstrate this, we printed UV-ICE circuits into polymer substrates with different mechanical properties to create resistive strain and force sensors; we printed flexible capacitive sensors with high sensitivity (2 kPa–1) and a wide range of measured pressures (from 5 Pa to 550 kPa) by creating a complex microstructure in the dielectric layer; we even realized ionic conductor-activated four-dimensional (4D) printing by printing a UV-ICE circuit into a shape memory polymer substrate. The proposed approach paves a new efficient way to realize multifunctional flexible devices and machines by bonding ICEs with other polymers in 3D forms.
  * <img src='/files/essay/he2022multimaterial.gif'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/he2022multimaterial.pdf)\] \[[Web](https://pubs.acs.org/doi/full/10.1021/acsami.2c18954)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/he2022multimaterial.txt)\]
  * Author: Xiangnan He, Jianxiang Cheng, Zhenqing Li, Haitao Ye, Xinfeng Wei, Honggeng Li, Rong Wang, Yuan-Fang Zhang, Hui Ying Yang, Chuanfei Guo, Qi Ge
  * Year: 2023
  * Journal: ACS Applied Materials & Interfaces
* Smart structures with embedded flexible sensors fabricated by fused deposition modeling-based multimaterial 3D printing
  * Keywords: 3D printing smart structures, flexible sensors
  * Summary: Smart structures have the advantages of high system integrity and diverse sensing capabilities. However, the labor-intensive and time-consuming fabrication process hinders the large-scale adoption of smart structures. Despite recent attempts to develop sensor-embedded structures using 3D printing technologies, the reported smart structures generally suffer from the complex fabrication process, constrained part size, and limited sensing modality. Herein, we propose a workflow to design and fabricate novel smart structures via multimaterial fused deposition modeling (FDM)-based 3D printing. More specifically, conductive filaments with tailorable mechanical and electrical properties, e.g. piezoresistive effects, were developed. Additionally, the printing process was optimized for processing soft filaments with Young’s modulus around 2 MPa, resolving the issue of filament buckling. Furthermore, the potential applications of the proposed workflow were showcased using three design cases, i.e. biaxial strain sensor, smart tire, and cable-driven soft finger with multiple sensing capabilities. This workflow provides a cost-effective and rapid solution for developing novel smart structures with soft materials.
  * <img src='/files/essay/ren2022smart.jpg'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/ren2022smart.pdf)\] \[[Web](https://www.tandfonline.com/doi/full/10.1080/19475411.2022.2095454)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/ren2022smart.txt)\]
  * Author: Huilin Ren, Xiaodan Yang, Zhenhu Wang, Xuguang Xu,R ong Wang, Qi Ge, Yi Xiong
  * Year: 2022
  * Journal: International Journal of Smart and Nano Materials
* Multi-jet ice 3D printing
  * Keywords: 3D printer development, Cryogenic additive manufacturing, Ice 3D printing, Multi-jet 3D printing
  * Summary: Multi-jet deposition of the materials is a matured technology used for graphic printing and 3 D printing for a wide range of materials. The multi-jet technology is fine-tuned for liquids with a specific range of viscosity and surface tension. However, the use of multi-jet for low viscosity fluids like water is not very popular. This paper aims to demonstrate the technique, particularly for the water-ice 3 D printing. 3 D printed ice parts can be used as patterns for investment casting, templates for microfluidic channel fabrication, support material for polymer 3 D printing, etc. Multi-jet ice 3 D printing is a novel technique for producing ice parts by selective deposition and freezing water layers. The paper confers the design, embodiment and integration of various subsystems of multi-jet ice 3 D printer. The outcomes of the machine trials are reported as case studies with elaborate details. The prismatic geometries are realized by ice 3 D printing. The accuracy of 0.1 mm is found in the build direction. The part height tends to increase due to volumetric expansion during the phase change. The present paper gives a novel architecture of the ice 3 D printer that produces the ice parts with good accuracy. The potential applications of the process are deliberated in this paper.
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/kamble2022multi.pdf)\] \[[Web](https://www.emerald.com/insight/content/doi/10.1108/RPJ-03-2021-0065/full/html?casa_token=Slpi9Fgpm30AAAAA:nz8355nNcRZBHklO726J5s4S9XQaPIXTKyde5ZWGjdhJVev6zNPWHdvX1_qz53bfUCu5h-VmLJ8oS8yRnBkM32ElgE3pTIpuOlXgXwpLMnWueUkbUGeUbQ)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/kamble2022multi.txt)\]
  * Author: Pushkar Prakash Kamble, Subodh Chavan, Rajendra Hodgir, Gopal Gote, K.P. Karunakaran
  * Year: 2022
  * Journal: Rapid Prototyping Journal
* Water/ice as sprayable sacrificial materials in low-temperature 3D printing for biomedical applications
  * Keywords: 3D printing, Additive manufacturing, Tissue engineering, Sacrificial process, Support structure, Low-temperature deposition
  * Summary: A spray-valve-based sacrificial process and “water/ice” support materials were developed in this study, based on the phase changes of water and low-temperature AM technology. The process entails utilizing the three phases of water to effectively and rapidly generate support structures that can also be efficiently removed completely after fabrication. Two scaffold materials, polycaprolactone-based waterborne polyurethane and chitosan, were tested. Moreover, ethanol or ethylene glycol was employed as the additive of pure water for facilitating the removal of excess support materials during fabrication. To verify the proposed sacrificial process, three complex and large scaffolds were fabricated, including a Y-type tubular scaffold, a square scaffold with dual inverted T-type inner channels, and a three-layer tapered scaffold. Results revealed that the properties of the printed complex structures were properly maintained and the adhesion between the layers was firm enough to resist elastic deformation.
  * <img src='/files/essay/liao2018water.jpg'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/liao2018water.pdf)\] \[[Web](https://www.sciencedirect.com/science/article/pii/S0264127518307640)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/liao2018water.txt)\]
  * Author: Chao-Yaug Liao, Wei-Jen Wu, Cheng-Tien Hsieh, Hung-Ching Yang, Ching-Shiow Tseng, Shan-hui Hsu
  * Year: 2018
  * Journal: Materials & Design
* Freeform 3D Ice Printing (3D-ICE) at the Micro Scale
  * Keywords: 3D Ice Printing
  * Summary: Water is one of the most important elements for life on earth. Water's rapid phase-change ability along with its environmental and biological compatibility also makes it a unique structural material for 3D printing of ice structures reproducibly and accurately. This work introduces the freeform 3D ice printing (3D-ICE) process for high-speed and reproducible fabrication of ice structures with micro-scale resolution. Drop-on-demand deposition of water onto a −35 °C platform rapidly transforms water into ice. The dimension and geometry of the structures are critically controlled by droplet ejection frequency modulation and stage motions. The freeform approach obviates layer-by-layer construction and support structures, even for overhang geometries. Complex and overhang geometries, branched hierarchical structures with smooth transitions, circular cross-sections, smooth surfaces, and micro-scale features (as small as 50 µm) are demonstrated. As a sample application, the ice templates are used as sacrificial geometries to produce resin parts with well-defined internal features. This approach could bring exciting opportunities for microfluidics, biomedical devices, soft electronics, and art.
  * <img src='/files/essay/garg2022freeform.jpg'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/garg2022freeform.pdf)\] \[[Web](https://onlinelibrary.wiley.com/doi/full/10.1002/advs.202201566)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/garg2022freeform.txt)\]
  * Author: Akash Garg, Saigopalakrishna S. Yerneni, Phil Campbell, Philip R. LeDuc, O. Burak Ozdoganlar
  * Year: 2022
  * Journal: Advanced Science
* Inkjet printing-based fabrication of microscale 3D ice structures
  * Keywords: microscale 3D printing ice structures
  * Summary: This study proposed a method for fabricating 3D microstructures of ice without a supporting material. The inkjet printing process was performed in a low humidity environment to precisely control the growth direction of the ice crystals. In the printing process, water droplets (volume = hundreds of picoliters) were deposited onto the previously formed ice structure, after which they immediately froze. Different 3D structures (maximum height = 2000 µm) could be formed by controlling the substrate temperature, ejection frequency and droplet size. The growth direction was dependent on the landing point of the droplet on the previously formed ice structure; thus, 3D structures could be created with high degrees of freedom.
  * <img src='/files/essay/zheng2020inkjet.webp'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/zheng2020inkjet.pdf)\] \[[Web](https://www.nature.com/articles/s41378-020-00199-x)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/zheng2020inkjet.txt)\]
  * Author: Fengyi Zheng, Zhongyan Wang, Jiasheng Huang, Zhihong Li
  * Year: 2020
  * Journal: Microsystems & Nanoengineering
* Ice lithography for 3D nanofabrication
  * Keywords:Nanotechnology, Nanofabrication, Electron-beam lithography, Ice lithography, 3D nanofabrication, Additive manufacturing, Organic ice
  * Summary: Nanotechnology and nanoscience are enabled by nanofabrication. Electron-beam lithography, which makes 2D patterns down to a few nanometers, is one of the fundamental pillars of nanofabrication. Recently, significant progress in 3D electron-beam-based nanofabrication has been made, such as the emerging ice lithography technology, in which ice thin-films are patterned by a focused electron-beam. Here, we review the history and progress of ice lithography, and focus on its applications in efficient 3D nanofabrication and additive manufacturing or nanoscale 3D printing. The finest linewidth made using frozen octane is below 5 nm, and nanostructures can be fabricated in selected areas on non-planar surfaces such as freely suspended nanotubes or nanowires. As developing custom instruments is required to advance this emerging technology, we discuss the evolution of ice lithography instruments and highlight major instrumentation advances. Finally, we present the perspectives of 3D printing of functional materials using organic ices. We believe that we barely scratched the surface of this new and exciting research area, and we hope that this review will stimulate cutting-edge and interdisciplinary research that exploits the undiscovered potentials of ice lithography for 3D photonics, electronics and 3D nanodevices for biology and medicine.
  * <img src='/files/essay/zhao2019ice.jpg'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/zhao2019ice.pdf)\] \[[Web](https://www.sciencedirect.com/science/article/pii/S209592731930324X)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/zhao2019ice.txt)\]
  * Author: Ding Zhao, Anpan Han, Min Qiu
  * Year: 2019
  * Journal: Science Bulletin
* Ice Lithography for Nanodevices
  * Keywords: Carbon nanotube e-beam lithography nanodevice field effect transistor
  * Summary: We report the successful application of a new approach, ice lithography (IL), to fabricate nanoscale devices. The entire IL process takes place inside a modified scanning electron microscope (SEM), where a vapor-deposited film of water ice serves as a resist for e-beam lithography, greatly simplifying and streamlining device fabrication. We show that labile nanostructures such as carbon nanotubes can be safely imaged in an SEM when coated in ice. The ice film is patterned at high e-beam intensity and serves as a mask for lift-off without the device degradation and contamination associated with e-beam imaging and polymer resist residues. We demonstrate the IL preparation of carbon nanotube field effect transistors with high-quality trans-conductance properties.
  * <img src='/files/essay/han2010ice.jpg'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/han2010ice.pdf)\] \[[Web](https://pubs.acs.org/doi/full/10.1021/nl1032815?casa_token=SOdhJB3z-mYAAAAA%3Ai6ESmLu1frRXJU0OpZ8T0Of-g1K5tzqY7fogXaKpoZl714o2Nk-jLMi_5av3KFROHekfOEwBCu7ial9naQ)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/han2010ice.txt)\]
  * Author: Anpan Han, Dimitar Vlassarev, Jenny Wang, Jene A. Golovchenko, Daniel Branton
  * Year: 2010
  * Journal: Nano letters
* Theoretical modeling of ice lithography on amorphous solid water
  * Keywords: ice lithography
  * Summary: Due to the perfection of the nanofabrication in nanotechnology and nanoscience, ice lithography (IL) by patterning ice thin-films with a focused electron beam, as a significant derivative technology of electron beam lithography (EBL), is attracting growing attention, evoked by its advantages over traditional EBL with respects of in situ-fabrication, high efficiency, high accuracy, limited proximity effect, three-dimensional (3D) profiling capability, etc. However, theoretical modeling of ice lithography for replicated profiles on the ice resist (amorphous solid water, ASW) has rarely been reported so far. As the result, the development of ice lithography still stays at the experimental stage. The shortage of modeling methods limits our insight into the ice lithography capability, as well as theoretical anticipations for future developments of this emerging technique. In this work, an e-beam induced etching ice model based on the Monte Carlo algorithm for point/line spread functions is established to calculate the replicated profiles of the resist by ice lithography. To testify the fidelity of the modeling method, systematic simulations of the ice lithography property under the processing parameters of the resist thickness, electron accelerating voltage and actual patterns are performed. Theoretical comparisons between the IL on ASW and the conventional EBL on polymethyl methacrylate (PMMA) show superior properties of IL over EBL in terms of the minimum feature size, the highest aspect ratio, 3D nanostructure/devices, etc. The success in developing a modeling method for ice lithography, as reported in this paper, offers a powerful tool in characterizing ice lithography up to the theoretical level and down to molecular scales.
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/liu2022theoretical.pdf)\] \[[Web](https://pubs.rsc.org/en/content/articlelanding/2022/nr/d2nr00594h/unauth)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/liu2022theoretical.txt)\]
  * Author: Tao Liu, Xujie Tong, Shuoqiu Tian, Yuying Xie, Mingsai Zhu, Bo Feng, Xiaohang Pan, Rui Zheng, Shan Wu, Ding Zhao, Yifang Chen, Bingrui Lu, Min Qiu
  * Year: 2022
  * Journal: Nanoscale
* 3D-Printed Phase-Change Artificial Muscles with Autonomous Vibration Control
  * Keywords: 3D printing, soft robot actuator, vibration control
  * Summary: Currently, additive manufacturing is utilized to fabricate many different actuators suited for soft robots. However, an effective controller paradigm is essential to benefit from the advantages of soft robots in terms of power consumption, production costs, weight, and safety while operating near living systems. In this work, an artificial muscle is additively manufactured with soft silicone elastomer material capable of demonstrating several levels of stiffness. The 3D-printed muscle is equipped with carbon fibers to receive a stimulus signal and develop a programmable joint that can present different stiffnesses. A nonlinear controller is developed to autonomously control the variable stiffness joint based on a reinforcement learning algorithm. The controller exhibits a slight increase in settling time; however, it demonstrates a decrease in fluctuation amplitude by 33% and a substantial reduction in power consumption by 41% in comparison to the optimized proportional integral derivative controller. At the same time, it is adaptable to and reliable in new conditions. The variable stiffness muscle is also used as a controllable mechanism to suppress the low frequency vibration. The study shows that the muscle can successfully attenuate the vibration autonomously when it is increased.
  * <img src='/files/essay/mohammadi20233d.jpg'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/mohammadi20233d.pdf)\] \[[Web](https://onlinelibrary.wiley.com/doi/full/10.1002/admt.202300199)\]
  * Author: Moslem Mohammadi, Abbas Z. Kouzani, Mahdi Bodaghi, Yong Xiang, Ali Zolfagharian
  * Year: 2023
  * Journal: Advanced Materials Technologies
* 3D printing of unsupported multi-scale and large-span ceramic via near-infrared assisted direct ink writing
  * Keywords: unsupported 3D Printing, Ceramic
  * Summary: In the three-dimensional printing process of ceramic with low-angle structures, additional supporting structures are usually employed to avoid collapse of overhanging parts. However, the extra supporting structures not only affect printing efficiency, but the problems caused by their removal are also a matter of concern. Herein, we present a ceramic printing method, which can realize printing of unsupported multi-scale and large-span ceramics through the combination of direct ink writing and near-infrared induced up-conversion particles-assisted photopolymerization. This printing technology enables in-situ curing of multi-scale filaments with diameters ranging from 410 µm to 3.50 mm, and ceramic structures of torsion spring, three-dimensional bending and cantilever beam were successfully constructed through unsupported printing. This method will bring more innovation to the unsupported 3D manufacturing of complex shape ceramics.
  * <img src='/files/essay/zhao20233d.webp'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/zhao20233d.pdf)\] \[[Web](https://www.nature.com/articles/s41467-023-38082-8)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/zhao20233d.txt)\]
  * Author: Yongqin Zhao, Junzhe Zhu, Wangyan He, Yu Liu, Xinxin Sang, Ren Liu
  * Year: 2023
  * Journal: Nature Communications
* Embedded 3D Printing of Architected Ceramics via Microwave-Activated Polymerization
  * Keywords: Ceramics 3D Printing
  * Summary: Light- and ink-based 3D printing methods have vastly expanded the design space and geometric complexity of architected ceramics. However, light-based methods are typically confined to a relatively narrow range of preceramic and particle-laden resins, while ink-based methods are limited in geometric complexity due to layerwise assembly. Here, embedded 3D printing is combined with microwave-activated curing to generate architected ceramics with spatially controlled composition in freeform shapes. Aqueous colloidal inks are printed within a support matrix, rapidly cured via microwave-activated polymerization, and subsequently dried and sintered into dense architectures composed of one or more oxide materials. This integrated manufacturing method opens new avenues for the design and fabrication of complex ceramic architectures with programmed composition, density, and form for myriad applications.
  * <img src='/files/essay/roman2023embedded.jpg'> <img src='/files/essay/roman2023embedded_1.jpg'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/roman2023embedded.pdf)\] \[[Web](https://onlinelibrary.wiley.com/doi/10.1002/adma.202209270)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/roman2023embedded.txt)\]
  * Author: Benito Román-Manso, Robert D. Weeks, Ryan L. Truby, Jennifer A. Lewis
  * Year: 2023
  * Journal: Advanced Materials
* Robot-assisted conformal additive manufacturing for continuous fibre-reinforced grid-stiffened shell structures
  * Keywords: Robotic system; grid stiffened shell structure; continuous fibre-reinforced polymer additive manufacturing; surface conformal toolpath
  * Summary: The advents in continuous fibre-reinforced polymer additive manufacturing (CFRP-AM) present unprecedented opportunities for the rapid development of next-generation high-performance composites with selectively and spatially distributed reinforcement. However, the widely adopted 3-degree-of-freedom motion configuration in current CFRP-AM systems hinders the exploration of composite structures with non-planar fibre layouts. This work presents a novel conformal CFRP-AM system to fabricate grid-stiffened shell structures leveraging its multi-DoF motion to pattern spatial features. The system integrates a 6-axis robot with an optimally designed co-extrusion module and operates through a design-to-manufacturing workflow. The proposed workflow includes three steps: system calibration, conformal toolpath generation, and process implementation. The conformal toolpath generation is a surface-mapping-based method that allows a simultaneous exploration of various geometric designs and their toolpaths. Experimental comparisons were made between parts fabricated by different processes, i.e., planar and conformal based, with different toolpaths, i.e., shells filled with zigzag and arc-offset patterns, and with various geometric designs, i.e., stiffener ribs with different crossline angles. The results manifest that the proposed system can significantly improve the compression strength and stiffness of grid-stiffened shell structures. Meanwhile, the additional design freedom on process and structure opens up a new possibility to customise their mechanical performance.
  * <img src='/files/essay/zhang2023robot.jpg'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/zhang2023robot.pdf)\] \[[Web](https://www.tandfonline.com/doi/full/10.1080/17452759.2023.2203695)\]
  * Author: Guoquan Zhang,Yaohui Wang,Ziwen Chen,Xuguang Xu,Ke Dong, Yi Xiong
  * Year: 2023
  * Journal: Virtual and Physical Prototyping
* 3D short fibre orientation for universal structures and geometries in material extrusion additive manufacturing
  * Keywords: 3D printing, Short fibre orientation, Fibre distribution, Fibre composites, Extruded filament geometry
  * Summary: Fibre orientation critically affects properties in material extrusion additive manufacturing. This study developed new understanding of 3D short-fibre orientation in representative structures intended to capture most features of typical tool-paths. The parametric representative structures included straight paths, curved paths, corners, a range of intersection types, and the start/end of paths. Fibre orientation tensors were used to quantify 3D fibre alignment along the direction of printed filaments (F-alignment), lateral to the in-plane filament direction (F-lat alignment), and normal to the print-platform (Z-alignment). Overall, fibres were highly aligned along the filament direction in both straight paths and curved paths but less aligned at corners and intersections. However, recovery of fibre orientation was found after corners and intersections. To assess fibre orientation uniformity throughout the layer thickness, specimens were sectioned normal to the print platform in seven planes throughout the thickness of a single extruded filament. High fibre orientation (F-alignment) was found intra-filament, but it gradually decreased from upper section to bottom section. Interlayer region showed both high F-alignment and Z-alignment. In addition, extruded-filament width and height critically affected fibre orientation: increasing extrusion width and layer height led to decreased F-alignment. Case studies showed the results are translatable to more complicated structures and a different polymer material. This study provides new understanding of 3D fibre orientation in additive manufacturing and will allow more informed design, analysis and optimisation of short-fibre-reinforced structures to improve performance.
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/yan20233d.pdf)\] \[[Web](https://www.sciencedirect.com/science/article/pii/S2214860423001483)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/yan20233d.txt)\]
  * Author: Jiongyi Yan, Emrah Demirci, Andrew Gleadall
  * Year: 2023
  * Journal: Additive Manufacturing
* Research and application of machine learning for additive manufacturing
  * Keywords: Review, Additive manufacturing, 3D Printing, Rapid prototyping, Machine learning, Deep learning, Digital manufacturing, Intelligent manufacturing
  * Summary: Additive manufacturing (AM) is poised to bring a revolution due to its unique production paradigm. It offers the prospect of mass customization, flexible production, on-demand and decentralized manufacturing. However, a number of challenges stem from not only the complexity of manufacturing systems but the demand for increasingly complex and high-quality products, in terms of design principles, standardization and quality control. These challenges build up barriers to the widespread adoption of AM in the industry and the in-depth research of AM in academia. To tackle the challenges, machine learning (ML) technologies rise to play a critical role as they are able to provide effective ways to quality control, process optimization, modelling of complex systems, and energy management. Hence, this paper employs a systematic literature review method as it is a defined and methodical way of identifying, assessing, and analysing published literature. Then, a keyword co-occurrence and cluster analysis are employed for analysing relevant literature. Several aspects of AM, including Design for AM (DfAM), material analytics, in situ monitoring and defect detection, property prediction and sustainability, have been clustered and summarized to present state-of-the-art research in the scope of ML for AM. Finally, the challenges and opportunities of ML for AM are uncovered and discussed.
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/qin2022research.pdf)\] \[[Web](https://doi.org/10.1016/j.addma.2022.102691)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/qin2022research.txt)\]
  * Author: Jian Qin, Fu Hu, Ying Liu, Paul Witherell, Charlie C.L. Wang, David W. Rosen, Timothy W. Simpson, Yan Lu, Qian Tang
  * Year: 2022
  * Journal: Additive Manufacturing
* Light-based vat-polymerization bioprinting
  * Keywords: vat-polymerization, review
  * Summary: Light-based vat-polymerization bioprinting enables computer-aided patterning of 3D cell-laden structures in a point-by-point, layer-by-layer or volumetric manner, using vat (vats) filled with photoactivatable bioresin (bioresins). This collection of technologies — divided by their modes of operation into stereolithography, digital light processing and volumetric additive manufacturing — has been extensively developed over the past few decades, leading to broad applications in biomedicine. In this Primer, we illustrate the methodology of light-based vat-polymerization 3D bioprinting from the perspectives of hardware, software and bioresin selections. We follow with discussions on methodological variations of these technologies, including their latest advancements, as well as elaborating on key assessments utilized towards ensuring qualities of the bioprinting procedures and products. We conclude by providing insights into future directions of light-based vat-polymerization methods.
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/levato2023light.pdf)\] \[[Web](https://doi.org/10.1038/s43586-023-00231-0)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/levato2023light.txt)\]
  * Author: Riccardo Levato, Oksana Dudaryeva, Carlos Ezio Garciamendez-Mijares, Bruce E. Kirkpatrick, Riccardo Rizzo, Jacob Schimelman, Kristi S. Anseth, Shaochen Chen, Marcy Zenobi-Wong, Yu Shrike Zhang
  * Year: 2023
  * Journal: Nature Reviews Methods Primers



<!-- *
  * Keywords:
  * Summary:
  * <img src='/files/essay/.'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/.pdf)\] \[[Web]()\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/.txt)\]
  * Author:
  * Year:
  * Journal:  -->

Soft Robotics
======
* Origami-based integration of robots that sense, decide, and respond
  * Keywords: origami robots by embedding sensing, computing, and actuating in compliant, conductive materials, without requiring semiconductor-based electronics
  * Summary: Origami-inspired engineering has enabled intelligent materials and structures to process and react to environmental stimuli. However, it is challenging to achieve complete sense-decide-act loops in origami materials for autonomous interaction with environments, mainly due to the lack of information processing units that can interface with sensing and actuation. Here, we introduce an integrated origami-based process to create autonomous robots by embedding sensing, computing, and actuating in compliant, conductive materials. By combining flexible bistable mechanisms and conductive thermal artificial muscles, we realize origami multiplexed switches and configure them to generate digital logic gates, memory bits, and thus integrated autonomous origami robots. We demonstrate with a flytrap-inspired robot that captures ‘living prey’, an untethered crawler that avoids obstacles, and a wheeled vehicle that locomotes with reprogrammable trajectories. Our method provides routes to achieve autonomy for origami robots through tight functional integration in compliant, conductive materials.
  * <img src='/files/essay/yan2023origami.webp'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/yan2023origami.pdf)\] \[[Web](https://www.nature.com/articles/s41467-023-37158-9#citeas)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/yan2023origami.txt)\]
  * Author: Wenzhong Yan, Shuguang Li, Mauricio Deguchi, Zhaoliang Zheng, Daniela Rus, Ankur Mehta
  * Year: 2023
  * Journal: Nature Communications
* On the Elastic Stability of Folded Rings in Circular and Straight States
  * Keywords: Folded Structures
  * Summary: Single-loop elastic rings can be folded into multi-loop equilibrium configurations. In this paper, the stability of several such multi-loop states which are either circular or straight are investigated analytically and illustrated by experimental demonstrations. The analysis ascertains stability by exploring variations of the elastic energy of the rings for admissible deformations in the vicinity of the equilibrium state. The approach employed is the conventional stability analysis for elastic conservative systems which differs from most of the analyses that have been published on this class of problems, as will be illustrated by reproducing and elaborating on several problems in the literature. In addition to providing solutions to two basic problems, the paper analyes and demonstrates the stability of six-sided rings that fold into straight configurations.
  * <img src='/files/essay/leanza2023elastic.jpg'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/leanza2023elastic.pdf)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/leanza2023elastic.txt)\]
  * Author: Sophie Leanza, Ruike Renee Zhao, John W. Hutchinson
  * Year: 2023
* Soft and lightweight fabric enables powerful and high-range pneumatic actuation
  * Keywords: a series of pneumatic actuators based on soft but less stretchable fabric, grasping force of over 150 N and a grasping range from 70 to 350 millimeters
  * Summary: Soft structures and actuation allow robots, conventionally consisting of rigid components, to perform more compliant, adaptive interactions similar to living creatures. Although numerous functions of these types of actuators have been demonstrated in the literature, their hyperelastic designs generally suffer from limited workspaces and load-carrying capabilities primarily due to their structural stretchability factor. Here, we describe a series of pneumatic actuators based on soft but less stretchable fabric that can simultaneously perform tunable workspace and bear a high payload. The motion mode of the actuator is programmable, combinable, and predictable and is informed by rapid response to low input pressure. A robotic gripper using three fabric actuators is also presented. The gripper demonstrates a grasping force of over 150 N and a grasping range from 70 to 350 millimeters. The design concept and comprehensive guidelines presented would provide design and analysis foundations for applying less stretchable yet soft materials in soft robots to further enhance their practicality.
  * <img src='/files/essay/zhang2023soft.jpg'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/zhang2023soft.pdf)\] \[[Web](https://www.science.org/doi/10.1126/sciadv.adg1203)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/zhang2023soft.txt)\]
  * Author: Zhuang Zhang, Yongzhou Long, Genliang Chen, Qichen Wu, Hao Wang, Hanqing Jiang
  * Year: 2023
  * Journal: Science Advances
* Predictive Learning of Error Recovery with a Sensorized Passivity-Based Soft Anthropomorphic Hand
  * Keywords: Soft robot sensing
  * Summary: Manipulation strategies based on the passive dynamics of soft-bodied interactions provide robust performances with limited sensory information. They utilize the kinematic structure and passive dynamics of the body to adapt to objects of varying shapes and properties. However, these soft passive interactions make the state of the robotic device influenced by the environment, making control generation and state estimation difficult. This work presents a closed-loop framework for dynamic interaction-based grasping that relies on two novelties: 1) a wrist-driven passive soft anthropomorphic hand that can generate robust grasp strategies using one-step kinaesthetic teaching and 2) a learning-based perception system that uses temporal data from sparse tactile sensors to predict and adapt to failures before it happens. With the anthropomorphic soft design and wrist-driven control, it is shown that controllers can be generated robust to novel objects and location uncertainty. With the learning-based high-level perception system and 32 sensing receptors, it is shown that failures can be predicted in advance, further improving the robustness of the entire system by more than doubling the grasping success rate. From over 1000 real-world grasping trials, both the control and perception framework are also seen to be transferable to novel objects and conditions. An interactive preprint version of the article can be found here:
  * <img src='/files/essay/gilday2023predictive.jpg'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/gilday2023predictive.pdf)\] \[[Web](https://onlinelibrary.wiley.com/doi/10.1002/aisy.202200390)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/gilday2023predictive.txt)\]
  * Author: Kieran Gilday, Thomas George-Thuruthel, Fumiya Iida
  * Year: 2023
  * Journal: Advanced Intelligent Systems
* Self-folding soft-robotic chains with reconfigurable shapes and functionalities
  * Keywords: Magnetic continuum soft robot
  * Summary: Magnetic continuum soft robots can actively steer their tip under an external magnetic field, enabling them to effectively navigate in complex in vivo environments and perform minimally invasive interventions. However, the geometries and functionalities of these robotic tools are limited by the inner diameter of the supporting catheter as well as the natural orifices and access ports of the human body. Here, we present a class of magnetic soft-robotic chains (MaSoChains) that can self-fold into large assemblies with stable configurations using a combination of elastic and magnetic energies. By pushing and pulling the MaSoChain relative to its catheter sheath, repeated assembly and disassembly with programmable shapes and functions are achieved. MaSoChains are compatible with state-of-the-art magnetic navigation technologies and provide many desirable features and functions that are difficult to realize through existing surgical tools. This strategy can be further customized and implemented for a wide spectrum of tools for minimally invasive interventions.
  * <img src='/files/essay/gu2023self.webp'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/gu2023self.pdf)\] \[[Video](https://ethz.ch/en/news-and-events/eth-news/news/2023/04/how-to-make-self-folding-surgical-tools.html)\] \[[Web](https://www.nature.com/articles/s41467-023-36819-z)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/gu2023self.txt)\]
  * Author: Hongri Gu, Marino Möckli, Claas Ehmke, Minsoo Kim, Matthias Wieland, Simon Moser, Clemens Bechinger, Quentin Boehler, Bradley J. Nelson
  * Year: 2023
  * Journal: Nature Communications
* A versatile jellyfish-like robotic platform for effective underwater propulsion and manipulation
  * Keywords: jellyfish-like robotic fish, high against-gravity speed, low input power
  * Summary: Underwater devices are critical for environmental applications. However, existing prototypes typically use bulky, noisy actuators and limited configurations. Consequently, they struggle to ensure noise-free and gentle interactions with underwater species when realizing practical functions. Therefore, we developed a jellyfish-like robotic platform enabled by a synergy of electrohydraulic actuators and a hybrid structure of rigid and soft components. Our 16-cm-diameter noise-free prototype could control the fluid flow to propel while manipulating objects to be kept beneath its body without physical contact, thereby enabling safer interactions. Its against-gravity speed was up to 6.1 cm/s, substantially quicker than other examples in literature, while only requiring a low input power of around 100 mW. Moreover, using the platform, we demonstrated contact-based object manipulation, fluidic mixing, shape adaptation, steering, wireless swimming, and cooperation of two to three robots. This study introduces a versatile jellyfish-like robotic platform with a wide range of functions for diverse applications.
  * <img src='/files/essay/wang2023versatile.jpg'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/wang2023versatile.pdf)\] \[[Web](https://www.science.org/doi/10.1126/sciadv.adg0292)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/wang2023versatile.txt)\]
  * Author: Tianlu Wang, Hyeong-Joon Joo, Shanyuan Song, Wenqi Hu, Christoph Keplinger, Metin Sitti
  * Year: 2023
  * Journal: Science Advances
* Kirigami-Inspired 3D Printable Soft Pneumatic Actuators with Multiple Deformation Modes for Soft Robotic Applications
  * Keywords: 3D printing, soft robots, Kirigami
  * Summary: In this article, a new soft pneumatic actuator (SPA) is proposed taking inspiration from Kirigami. Kirigami-inspired cuts are applied to the actuator design, which enables the SPA to be equipped with multiple deformation modes. The proposed Kirigami-inspired soft pneumatic actuator (KiriSPA) is capable of producing bending motion, stretching motion, contraction motion, combined motion of bending and stretching, and combined motion of bending and contraction. The KiriSPA can be directly manufactured using 3D printers based on the fused deposition modeling technology. Finite element method is used to analyze and predict the deformation modes of the KiriSPA. We also investigated the step response, creep, hysteresis, actuation speed, stroke, workspace, stiffness, power density, and blocked force of the KiriSPA. Moreover, we demonstrated that KiriSPAs can be combined to expand the capabilities of various soft robotic systems including the soft robotic gripper for delicate object manipulation, the soft planar robotic manipulator for picking objects in the confined environment, the quadrupedal soft crawling robot, and the soft robot with the flipping locomotion.
  * <img src='/files/essay/guo2023kirigami.jpg'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/guo2023kirigami.pdf)\] \[[Web](https://www.liebertpub.com/doi/10.1089/soro.2021.0199)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/guo2023kirigami.txt)\]
  * Author: Jin Guo, Zeyu Li, Jin-Huat Low, Qianqian Han, Chao-Yu Chen, Jun Liu, Zhuangjian Liu, Chen-Hua Yeow
  * Year: 2023
  * Journal: Soft Robotics
* Design Optimization of Soft Robots: A Review of the State of the Art
  * Keywords:
  * Summary: Robotics has undergone a profound revolution in the past 50 years, moving from the laboratory and research institute to the factory and home. Kinematics and dynamics theories have been developed as the foundation for robot design and control, based on the conventional definition of robots: a kinematic chain of rigid links.
  * <img src='/files/essay/chen2020design.gif'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/chen2020design.pdf)\] \[[Web](https://ieeexplore.ieee.org/document/9237112)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/chen2020design.txt)\]
  * Author: Feifei Chen, Michael Yu Wang
  * Year: 2020
  * Journal: IEEE Robotics & Automation Magazine
* Advanced Soft Robotic System for In Situ 3D Bioprinting and Endoscopic Surgery
  * Keywords: Soft Robotic, 3D Bioprinting, Surgery
  * Summary: Three-dimensional (3D) bioprinting technology offers great potential in the treatment of tissue and organ damage. Conventional approaches generally rely on a large form factor desktop bioprinter to create in vitro 3D living constructs before introducing them into the patient's body, which poses several drawbacks such as surface mismatches, structure damage, and high contamination along with tissue injury due to transport and large open-field surgery. In situ bioprinting inside a living body is a potentially transformational solution as the body serves as an excellent bioreactor. This work introduces a multifunctional and flexible in situ 3D bioprinter (F3DB), which features a high degree of freedom soft printing head integrated into a flexible robotic arm to deliver multilayered biomaterials to internal organs/tissues. The device has a master-slave architecture and is operated by a kinematic inversion model and learning-based controllers. The 3D printing capabilities with different patterns, surfaces, and on a colon phantom are also tested with different composite hydrogels and biomaterials. The F3DB capability to perform endoscopic surgery is further demonstrated with fresh porcine tissue. The new system is expected to bridge a gap in the field of in situ bioprinting and support the future development of advanced endoscopic surgical robots.
  * <img src='/files/essay/thai2023advanced.jpg'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/thai2023advanced.pdf)\] \[[Web](https://onlinelibrary.wiley.com/doi/10.1002/advs.202205656)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/thai2023advanced.txt)\]
  * Author: Mai Thanh Thai, Phuoc Thien Phan, Hien Anh Tran, Chi Cong Nguyen, Trung Thien Hoang, James Davies, Jelena Rnjak-Kovacina, Hoang-Phuong Phan, Nigel Hamilton Lovell, Thanh Nho Do
  * Year: 2023
  * Journal: Advanced Science
* Liquid metal droplets bouncing higher on thicker water layer
  * Keywords: Liquid metal
  * Summary: Liquid metal (LM) has gained increasing attention for a wide range of applications, such as flexible electronics, soft robots, and chip cooling devices, owing to its low melting temperature, good flexibility, and high electrical and thermal conductivity. In ambient conditions, LM is susceptible to the coverage of a thin oxide layer, resulting in unwanted adhesion with underlying substrates that undercuts its originally high mobility. Here, we discover an unusual phenomenon characterized by the complete rebound of LM droplets from the water layer with negligible adhesion. More counterintuitively, the restitution coefficient, defined as the ratio between the droplet velocities after and before impact, increases with water layer thickness. We reveal that the complete rebound of LM droplets originates from the trapping of a thinly low-viscosity water lubrication film that prevents droplet-solid contact with low viscous dissipation, and the restitution coefficient is modulated by the negative capillary pressure in the lubrication film as a result of the spontaneous spreading of water on the LM droplet. Our findings advance the fundamental understanding of complex fluids’ droplet dynamics and provide insights for fluid control.
  * <img src='/files/essay/dai2023liquid.webp'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/dai2023liquid.pdf)\] \[[Web](https://doi.org/10.1038/s41467-023-39348-x)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/dai2023liquid.txt)\]
  * Author: Yuhang Dai, Minfei Li, Bingqiang Ji, Xiong Wang, Siyan Yang, Peng Yu, Steven Wang, Chonglei Hao, Zuankai Wang
  * Year: 2023
  * Journal: Nature Communications
* Desktop fabrication of monolithic soft robotic devices with embedded fluidic control circuits
  * Keywords: 3D printing pneumatic soft robots
  * Summary: Most soft robots are pneumatically actuated and fabricated by molding and assembling processes that typically require many manual operations and limit complexity. Furthermore, complex control components (for example, electronic pumps and microcontrollers) must be added to achieve even simple functions. Desktop fused filament fabrication (FFF) three-dimensional printing provides an accessible alternative with less manual work and the capability of generating more complex structures. However, because of material and process limitations, FFF-printed soft robots often have a high effective stiffness and contain a large number of leaks, limiting their applications. We present an approach for the design and fabrication of soft, airtight pneumatic robotic devices using FFF to simultaneously print actuators with embedded fluidic control components. We demonstrated this approach by printing actuators an order of magnitude softer than those previously fabricated using FFF and capable of bending to form a complete circle. Similarly, we printed pneumatic valves that control a high-pressure airflow with low control pressure. Combining the actuators and valves, we demonstrated a monolithically printed electronics-free autonomous gripper. When connected to a constant supply of air pressure, the gripper autonomously detected and gripped an object and released the object when it detected a force due to the weight of the object acting perpendicular to the gripper. The entire fabrication process of the gripper required no posttreatment, postassembly, or repair of manufacturing defects, making this approach highly repeatable and accessible. Our proposed approach represents a step toward complex, customized robotic systems and components created at distributed fabricating facilities.
  * <img src='/files/essay/zhai2023desktop.jpg'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/zhai2023desktop.pdf)\] \[[Web](https://doi.org/10.1126/scirobotics.adg3792)\]
  * Author: Yichen Zhai, Albert De Boer, Jiayao Yan, Benjamin Shih, Martin Faber, Joshua Speros, Rohini Gupta, Michael T. Tolley
  * Year: 2023
  * Journal: Science Robotics
* Caterpillar-inspired soft crawling robot with distributed programmable thermal actuation
  * Keywords: Caterpillar, soft robots
  * Summary: Many inspirations for soft robotics are from the natural world, such as octopuses, snakes, and caterpillars. Here, we report a caterpillar-inspired, energy-efficient crawling robot with multiple crawling modes, enabled by joule heating of a patterned soft heater consisting of silver nanowire networks in a liquid crystal elastomer (LCE)–based thermal bimorph actuator. With patterned and distributed heaters and programmable heating, different temperature and hence curvature distribution along the body of the robot are achieved, enabling bidirectional locomotion as a result of the friction competition between the front and rear end with the ground. The thermal bimorph behavior is studied to predict and optimize the local curvature of the robot under thermal stimuli. The bidirectional actuation modes with the crawling speeds are investigated. The capability of passing through obstacles with limited spacing are demonstrated. The strategy of distributed and programmable heating and actuation with thermal responsive materials offers unprecedented capabilities for smart and multifunctional soft robots.
  * <img src='/files/essay/wu2023caterpillar.jpg'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/wu2023caterpillar.pdf)\] \[[Web](https://doi.org/10.1126/sciadv.adf8014)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/wu2023caterpillar.txt)\]
  * Author: Shuang Wu, Yaoye Hong, Yao Zhao, Jie Yin, Yong Zhu
  * Year: 2023
  * Journal: Science Advances


<!-- *
  * Keywords:
  * Summary:
  * <img src='/files/essay/.'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/.pdf)\] \[[Web](https://doi.org/)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/.txt)\]
  * Author:
  * Year:
  * Journal:  -->

Smart Materials & Structures
======
* A stretchable wireless wearable bioelectronic system for multiplexed monitoring and combination treatment of infected chronic wounds
  * Keywords: wearable bioelectronic system
  * Summary: Chronic nonhealing wounds are one of the major and rapidly growing clinical complications all over the world. Current therapies frequently require emergent surgical interventions, while abuse and misapplication of therapeutic drugs often lead to an increased morbidity and mortality rate. Here, we introduce a wearable bioelectronic system that wirelessly and continuously monitors the physiological conditions of the wound bed via a custom-developed multiplexed multimodal electrochemical biosensor array and performs noninvasive combination therapy through controlled anti-inflammatory antimicrobial treatment and electrically stimulated tissue regeneration. The wearable patch is fully biocompatible, mechanically flexible, stretchable, and can conformally adhere to the skin wound throughout the entire healing process. Real-time metabolic and inflammatory monitoring in a series of preclinical in vivo experiments showed high accuracy and electrochemical stability of the wearable patch for multiplexed spatial and temporal wound biomarker analysis. The combination therapy enabled substantially accelerated cutaneous chronic wound healing in a rodent model.
  * <img src='/files/essay/shirzaei2023stretchable.jpg'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/shirzaei2023stretchable.pdf)\] \[[Web](https://www.science.org/doi/10.1126/sciadv.adf7388)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/shirzaei2023stretchable.txt)\]
  * Author: Ehsan Shirzaei Sani, Changhao Xu, Canran Wang, Yu Song, Jihong Min, Jiaobing Tu, Samuel A. Solomon, Jiahong Li, Jaminelli L. Banks, David G. Armstrong, Wei Gao
  * Year: 2023
  * Journal: Science Advances
* Investigating the effect of surface protrusions on galloping energy harvesting Editor’s Pick
  * Keywords: piezoelectric energy harvester
  * Summary: This Letter explores the potential effect of implementing different surface protrusions on galloping energy harvesters. Three types of protruded bluff bodies with rectangular, triangular, and elliptical metasurfaces are proposed, and four kinds of surface treatments are deployed to vary their protruded shape. Wind tunnel experiments reveal that adding the protrusions can obviously change the mode of oscillations, and only the backward protrusions can enhance the galloping response. Both the experiments and simulations show that elliptical surface protrusions have the greatest potential to enhance the galloping energy harvesting performance. Specifically, with a backward protruded length of 15 mm, the maximum output power in the experiments is measured to be 0.757 mW, which occurs at 5.1 m/s, and an optimal load resistance of 300 kΩ. In this case, the energy harvester outperforms its counterpart carrying a simple square prism by 157.48%.
  * <img src='/files/essay/xing2023investigating.jpg'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/xing2023investigating.pdf)\] \[[Web](https://pubs.aip.org/aip/apl/article/122/15/153902/2878606/Investigating-the-effect-of-surface-protrusions-on)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/xing2023investigating.txt)\]
  * Author: Juntong Xing, Masoud Rezaei, Huliang Dai, Wei-Hsin Liao
  * Year: 2023
  * Journal: Applied Physics Letters
* Tailorable activation of thermoresponsive composite structures incorporating wavy heaters via hybrid manufacturing
  * Keywords: Joule-heating, Tailorable activation of smart materials, Shape-memory polymers, Thermoresponsive composite structures
  * Summary: Herein, we present a concept of embedding a wavy heater into a thermoresponsive material matrix to form a composite structure with parametrically designed thermal activation behavior, through a facile manufacturing approach combining 3D-printing and laser-cutting. We develop a numerical model to predict the transient heat transfer for varying wavy shapes of heater, and experimentally validate the numerical results. The exploration of the design space using the numerical model shows a reduction of up to 82% in heating time using the wavy design compared with the flat design. Finally, we demonstrate experimentally the stiffness tuning in thermoresponsive composite structures. This work paves the way for large-scale thermoresponsive composite structures with applications in aerospace and architecture.
  * <img src='/files/essay/zhang2023tailorable.jpg'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/zhang2023tailorable.pdf)\] \[[Web](https://www.sciencedirect.com/science/article/pii/S2452213923000311)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/zhang2023tailorable.txt)\]
  * Author: Yuan-Fang Zhang, Honggeng Li, Chengyun Long, Yi Xiong, Qi Ge
  * Year: 2023
  * Journal: Composites Communications

Metamaterials
=====
* Mechanical metamaterials made of freestanding quasi-BCC nanolattices of gold and copper with ultra-high energy absorption capacity
  * Keywords: Mechanical metamaterials, ultra-high energy absorption capacity, nm-scale manufacturing & testing
  * Summary: Nanolattices exhibit attractive mechanical properties such as high strength, high specific strength, and high energy absorption. However, at present, such materials cannot achieve effective fusion of the above properties and scalable production, which hinders their applications in energy conversion and other fields. Herein, we report gold and copper quasi-body centered cubic (quasi-BCC) nanolattices with the diameter of the nanobeams as small as 34 nm. We show that the compressive yield strengths of quasi-BCC nanolattices even exceed those of their bulk counterparts, despite their relative densities below 0.5. Simultaneously, these quasi-BCC nanolattices exhibit ultrahigh energy absorption capacities, i.e., 100 ± 6 MJ m−3 for gold quasi-BCC nanolattice and 110 ± 10 MJ m−3 for copper quasi-BCC nanolattice. Finite element simulations and theoretical calculations reveal that the deformation of quasi-BCC nanolattice is dominated by nanobeam bending. And the anomalous energy absorption capacities substantially stem from the synergy of the naturally high mechanical strength and plasticity of metals, the size reduction-induced mechanical enhancement, and the quasi-BCC nanolattice architecture. Since the sample size can be scaled up to macroscale at high efficiency and affordable cost, the quasi-BCC nanolattices with ultrahigh energy absorption capacity reported in this work may find great potentials in heat transfer, electric conduction, catalysis applications.
  * <img src='/files/essay/cheng2023mechanical.webp'> <img src='/files/essay/cheng2023mechanical_2.webp'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/cheng2023mechanical.pdf)\] \[[Web](https://www.nature.com/articles/s41467-023-36965-4)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/cheng2023mechanical.txt)\]
  * Author: Hongwei Cheng, Xiaoxia Zhu, Xiaowei Cheng, Pengzhan Cai, Jie Liu, Huijun Yao, Ling Zhang, Jinglai Duan
  * Year: 2023
  * Journal: Nature Communications
* Multimaterial 3D printed self-locking thick-panel origami metamaterials
  * Keywords: 3D self-locking thick-panel origami structure
  * Summary: Thick-panel origami has shown great potential in engineering applications. However, the thick-panel origami created by current design methods cannot be readily adopted to structural applications due to the inefficient manufacturing methods. Here, we report a design and manufacturing strategy for creating thick-panel origami structures with excellent foldability and capability of withstanding cyclic loading. We directly print thick-panel origami through a single fused deposition modeling (FDM) multimaterial 3D printer following a wrapping-based fabrication strategy where the rigid panels are wrapped and connected by highly stretchable soft parts. Through stacking two thick-panel origami panels into a predetermined configuration, we develop a 3D self-locking thick-panel origami structure that deforms by following a push-to-pull mode enabling the origami structure to support a load over 11000 times of its own weight and sustain more than 100 cycles of 40% compressive strain. After optimizing geometric parameters through a self-built theoretical model, we demonstrate that the mechanical response of the self-locking thick-panel origami structure is highly programmable, and such multi-layer origami structure can have a substantially improved impact energy absorption for various structural applications.
  * <img src='/files/essay/ye2023multimaterial.webp'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/ye2023multimaterial.pdf)\] \[[Web](https://www.nature.com/articles/s41467-023-37343-w)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/ye2023multimaterial.txt)\]
  * Author: Haitao Ye, Qingjiang Liu, Jianxiang Cheng, Honggeng Li, Bingcong Jian, Rong Wang, Zechu Sun, Yang Lu, Qi Ge
  * Year: 2023
  * Journal: Nature Communications
* Encoding and Storage of Information in Mechanical Metamaterials
  * Keywords: Metamaterials, encoding and storage information
  * Summary: Information processing using material's own properties has gained increasing interest. Mechanical metamaterials, due to their diversity of deformation modes and wide design space, can be used to realize information processing, such as computing and storage. Here a mechanical metamaterial system is demonstrated for material-based encoding and storage of data through programmed reconfigurations of the metamaterial's structured building blocks. Sequential encoding and decoding are achieved in the three-dimensional (3D) printed pixelated mechanical metamaterial via kirigami-based “pixels” with programmable, temperature-dependent bistability. The mechanical metamaterial is demonstrated via a multistep deformation of encoding messages of texts and surfaces with arrays of binary data, and then decoding them by applying a predetermined stretching and heating regimen to sequentially retrieve layers of stored information and display them on its surface. This approach serves as a general framework to enable the encoding and storage of data with mechanical metamaterials.
  * <img src='/files/essay/meng2023encoding.jpg'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/meng2023encoding.pdf)\] \[[Web](https://onlinelibrary.wiley.com/doi/full/10.1002/advs.202301581)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/meng2023encoding.txt)\]
  * Author: Zhiqiang Meng, Hujie Yan, Mingchao Liu, Wenkai Qin, Guy M. Genin, Chang Qing Chen
  * Year: 2023
  * Journal: Advanced Science
* Design of a programmable particle filtering medium using a novel auxetic metamaterial
  * Keywords: Auxetic, Metamaterial, Additive manufacturing, Particle filter, Programmable, Genetic algorithm, Optimization
  * Summary: This manuscript aims to design and develop a 2D auxetic filtering medium with programmable geometric features specifically designed to vary under in-plane tensile strain. This feature empowers the filtering medium to control the particles separation. A novel design and optimisation algorithm developed in Matlab® determines the final optimized geometry of the filtering medium based on the desired particle size input. Upon thorough numerical investigation, an empirical relationship between the linear elastic in-plane tensile strain and aperture size of the proposed metamaterial is revealed. This empirical relation can be used in mechatronic and control systems to steer the proposed filtering medium. A prototype of such filtering medium capable of classification of particles of size 4mm to 4.5mm, when subjected to linear strain, is fabricated through Fused Deposition Modelling (FDM) process. The developed geometry configurations in this research are scalable, providing a potential cost-effective and efficient solution for industrial applications including reconfigurable filtration and segregation systems.
  * <img src='/files/essay/ali2023design.png'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/ali2023design.pdf)\] \[[Web](https://iopscience.iop.org/article/10.1088/1361-665X/acceea/meta)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/ali2023design.txt)\]
  * Author: Hafiz Muhammad Asad Ali, Meisam Abdi, Sayyed Abolfazl Zahedi, Yong Sun
  * Year: 2023
  * Journal: Smart Materials and Structures

<!-- *
  * Keywords:
  * Summary:
  * <img src='/files/essay/.'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/.pdf)\] \[[Web]()\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/.txt)\]
  * Author:
  * Year:
  * Journal:  -->


Topology Optimization
=====
* Interactive Structural Topology Optimization with Subjective Scoring and Drawing Systems
  * Keywords: Topology optimization, Subjective preferences, Scoring, Drawing, Structural design
  * Summary: Topology optimization techniques can create efficient and innovative structural designs by redistributing underutilized materials to the most-needed locations. These techniques are typically performed based purely on structural performance without considering factors like aesthetics and other design requirements. Hence, the obtained structural designs may not be suitable for specific practical applications. This study presents a new topology optimization method, SP-BESO, by considering the subjective preferences (SP) of the designer. Here, subjective scoring and drawing systems are introduced into the popular bi-directional evolutionary structural optimization (BESO) technique. The proposed SP-BESO method allows users to iteratively and interactively create topologically different and structurally efficient solutions by explicitly scoring and drawing their subjective preferences. Hence, users do not need to passively accept the optimization results. A user-friendly digital design tool, iBESO, is developed, which contains four optimizers to simultaneously perform the proposed SP-BESO method to assist in the design exploration task. A variety of 2D examples are tested using the iBESO software to demonstrate the effectiveness of the proposed SP-BESO method. It is found that the combination of parameters used in the scoring and drawing systems controls the formation of final structural topologies toward performance-driven or preference-driven designs. The utilization of the proposed SP-BESO method in potential practical applications is also demonstrated.
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/li2023interactive.pdf)\] \[[Video](https://www.linkedin.com/posts/yi-min-mike-xie-11015b118_topologyoptimization-preferences-structuraldesign-ugcPost-7052119601622519808-4UC7?utm_source=share&utm_medium=member_desktop)\] \[[Web](https://www.sciencedirect.com/science/article/pii/S0010448523000647?dgcid=author)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/li2023interactive.txt)\]
  * Author: Zhi Li, Ting-Uei Lee, Yi Min Xie
  * Year: 2023
  * Journal: Computer-Aided Design
* 3D-printed high-toughness composite structures by anisotropic topology optimization
  * Keywords: Polymer-matrix composites (PMCs), Carbon fibers, Anisotropy3D printing
  * Summary: The toughness of structures is essential to prevent catastrophic failure. This study introduced a design framework to improve the toughness of 3D-printed carbon fiber-reinforced composite structures by local latticing utilizing the intermediate material fraction obtained in the topology optimization. The framework was based on anisotropic topology optimization considering material fraction and material orientation. The optimized results were de-homogenized by the phase field-based technique to determine the 3D printing path. Experimental validations were carried out on a three-point bending beam problem. As a result, it was shown that the framework endowed toughness for the 3D-printed carbon fiber-reinforced composite structure.
  * <img src='/files/essay/ichihara20233d_1.jpg'> <img src='/files/essay/ichihara20233d_2.jpg'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/ichihara20233d.pdf)\] \[[Web](https://www.sciencedirect.com/science/article/pii/S1359836823000756)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/ichihara20233d.txt)\]
  * Author: Naruki Ichihara, Masahito Ueda
  * Year: 2023
  * Journal: Composites Part B: Engineering
* Double-Spirals Offer the Development of Preprogrammable Modular Metastructures
  * Keywords: Preprogrammable Modular Metastructures
  * Summary: Metamaterials with adjustable, sometimes unusual properties offer advantages over conventional materials with predefined mechanical properties in many technological applications. A group of metamaterials, called modular metamaterials or metastructures, are developed through the arrangement of multiple, mostly similar building blocks. These modular structures can be assembled using prefabricated modules and reconfigured to promote efficiency and functionality. Herein, a novel modular metastructure is developed by taking advantage of the high compliance of preprogrammable double-spirals. First, the mechanical behavior of a four-module metastructure under tension, compression, rotation, and sliding is simulated using the finite-element method. Then, 3D printing and mechanical testing are used to illustrate the tunable anisotropic and asymmetric behavior of the spiral-based metastructures in practice. The results show the simple reconfiguration of the presented metastructure toward the desired functions. The mechanical behavior of single double-spirals and the characteristics that can be achieved through their combinations make our modular metastructure suitable for various applications in robotics, aerospace, and medical engineering.
  * <img src='/files/essay/jafarpour2023double.jpg'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/jafarpour2023double.pdf)\] \[[Web](https://onlinelibrary.wiley.com/doi/full/10.1002/adem.202300102)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/jafarpour2023double.txt)\]
  * Author: Mohsen Jafarpour, Stanislav N. Gorb, Hamed Rajabi
  * Year: 2023
  * Journal: Advanced Engineering Materials
* Active Materials for Functional Origami
  * Keywords: Active Materials, Functional Origami, Review
  * Summary: In recent decades, origami has been explored to aid in the design of engineering structures. These structures span multiple scales and have been demonstrated to be used towards various areas such as aerospace, metamaterial, biomedical, robotics, and architectural applications. Conventionally, origami or deployable structures have been actuated by hands, motors, or pneumatic actuators, which can result in heavy or bulky structures. On the other hand, active materials, which reconfigure in response to external stimulus, eliminate the need for external mechanical loads and bulky actuation systems. Thus, in recent years, active materials incorporated with deployable structures have shown promise for remote actuation of light weight, programmable origami. In this review, active materials such as shape memory polymers and alloys, hydrogels, liquid crystal elastomers, magnetic soft materials, and covalent adaptable network polymers, their actuation mechanisms, as well as how they have been utilized for active origami and where these structures are applicable is discussed. Additionally, the state-of-the-art fabrication methods to construct active origami are highlighted. The existing structural modeling strategies for origami, the constitutive models used to describe active materials, and the largest challenges and future directions for active origami research are summarized.
  * <img src='/files/essay/leanza2023active.png'>
  * \[[PDF](http://Liuchao-JIN.github.io/files/essay/leanza2023active.pdf)\] \[[Web](https://onlinelibrary.wiley.com/doi/abs/10.1002/adma.202302066)\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/leanza2023active.txt)\]
  * Author: Sophie Leanza, Shuai Wu, Xiaohao Sun, H. Jerry Qi, Ruike Renee Zhao
  * Year: 2023
  * Journal: Advanced Materials


  <!-- *
    * Keywords:
    * Summary:
    * <img src='/files/essay/.'>
    * \[[PDF](http://Liuchao-JIN.github.io/files/essay/.pdf)\] \[[Web]()\] \[[BibTeX](http://Liuchao-JIN.github.io/files/essay/.txt)\]
    * Author:
    * Year:
    * Journal:  -->
