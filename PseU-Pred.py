import streamlit as st
from streamlit_option_menu import option_menu
from PIL import Image

selected2 = option_menu(None, ["Home", "Predictor", "Dataset", "Citations"],
                        icons=['house', 'search', 'list-task', 'book'],
                        menu_icon="cast", default_index=0, orientation="horizontal")
# selected2

if selected2 == "Home":
    st.title('PseU-Pred')
    #st.subheader("The DHU-Pred is a web-server for the prediction of Dihydrouridine in transfer RNA (tRNA) "
     #            "modifications.")
    #st.write("The current study focused on the most prevalent uridine modification, named pseudouridine (ψ)."
     #        "Uridine is widely produced in nature as uridine monophosphate (uridylate) by de novo synthesis by the decarboxylation of orotidylate, which is catalysed by orotidylate decarboxylase."
      #       "The current research work focused on the identification of ψ sites by extracting features from the RNA sequences using novel feature development methodology. "
       #      "The benchmark dataset obtained from chen et al (Chen et al., 2016) were used for training as well as model testing and validation."
        #     "The dataset composed of RNA sequences belonging to three species homosapiens, saccharomyces cerevisiae and mus musculus designated as HS_990, SC_628 and MM_944 respectively."
         #    )
    image = Image.open('Flowchart.png')
    st.image(image, width=400)

elif selected2 == "Predictor":
    #st.subheader("Predictor Page")
    import predictor

    exec(open('predictor.py').read())


elif selected2 == "Dataset":
    #st.subheader("Data Set")
    st.info("Homosapiens (HSP_990) FASTA file")
    with open("HSP_990.fas", "rb") as file:
        btn = st.download_button(
            label="Download file",
            data=file,
            file_name="HSP_990.fas",
            mime=""
        )
    st.info("Saccharomyces cerevisiae (SCV_628) FASTA file")
    with open("SCV_628.fas", "rb") as file:
        btn = st.download_button(
            label="Download file",
            data=file,
            file_name="SCV_628.fas",
            mime=""
        )
    st.info("Mus Musculus (MMC_944) FASTA file")
    with open("MMC_944.fas", "rb") as file:
        btn = st.download_button(
            label="Download file",
            data=file,
            file_name="MMC_944.fas",
            mime=""
        )
    st.info("Indpendent Homosapiens (HSP_200) FASTA file")
    with open("HSP_Independent.fas", "rb") as file:
        btn = st.download_button(
            label="Download file",
            data=file,
            file_name="HSP_Independent.fas",
            mime=""
        )
    st.info("Indpendent Saccahromyces cerevisiae (SCV_200) FASTA file")
    with open("SCV_Independent.fas", "rb") as file:
        btn = st.download_button(
            label="Download file",
            data=file,
            file_name="SCV_Independent.fas",
            mime=""
        )
