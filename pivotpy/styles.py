light_style = """<style>
               .widget-text input { 
                   background-color:white !important;
                   border-radius:20px !important; 
                   padding: 0px 10px 0px 10px !important;
                   border: 1px solid whitesmoke !important;
                   color: gray !important;
                   }
               .widget-text input:focus {
                   border: 1px solid skyblue !important;
                   }
               .widget-text input:hover { 
                   border: 1px solid skyblue !important;
                   }
               .widget-dropdown > select {
                   background-color: transparent !important;
                   border:none !important;
                   border-bottom: 1px solid skyblue !important;
                   box-shadow: inset 0 -20px 10px -20px skyblue;
               }
               .widget-dropdown > select:hover {
                   background-color: white !important;
               }
               .widget-dropdown > select > option {
                   color: gray !important;
                   background-color: #eaf0f0 !important;
               }
               .widget-dropdown > select > option:focus {
                   background-color: red !important;
               }
               .widget-box { 
                   background-color:#eaf0f0 !important;
                   border-radius:5px !important; 
                   padding:10px !important;
                   border: 1px solid whitesmoke !important;
                   box-shadow: 1px 1px 1px 1px gray !important;
                }
                .borderless {
                border: 1px solid transparent !important;
                box-shadow: none !important;
                border-radius: 0px !important;
                margin:4px !important;
                }
                .marginless {
                    margin: 0px !important;
                    border-radius: 0px !important;
                }
                .widget-tab {
                   background-color:#eaf0f0 !important;  
                   border: 1px solid whitesmoke !important;
                   box-shadow: 1px 1px 1px 1px gray !important;
                }
            </style>"""

dark_style = """<style>
               .widget-text input { 
                   background-color:black !important;
                   border-radius:20px !important; 
                   padding: 0px 10px 0px 10px !important;
                   border: 1px solid skyblue !important;
                   color: #abb2bf !important;
                   }
               .widget-text input:focus {
                   border: 1px solid skyblue !important;
                   }
               .widget-text input:hover { 
                   border: 1px solid skyblue !important;
                   }
               .widget-dropdown > select {
                   background-color: transparent !important;
                   border:none !important;
                   border-bottom: 1px solid skyblue !important;
                   box-shadow: inset 0 -20px 10px -20px skyblue;
                   color: white !important;
               }
               .widget-dropdown > select:hover {
                   background-color: black !important;
               }
               .widget-dropdown > select > option {
                   color: whitesmoke !important;
                   background-color: black !important;
               }
               .widget-dropdown > select > option:focus {
                   background-color: red !important;
               }
               .widget-label {
                   color: white !important;
               }
               .widget-html {
                   color: white !important;
               }
               .widget-box { 
                   background-color:#282c34 !important;
                   border-radius:5px !important; 
                   padding:10px !important;
                   border: 1px solid black !important;
                   box-shadow: 1px 1px 1px 1px black !important;
                }
                .borderless {
                border: 1px solid transparent !important;
                box-shadow: none !important;
                border-radius: 4px !important;
                margin:4px !important;
                }
                .marginless {
                    margin: 0px !important;
                    border-radius: 0px !important;
                }
                .widget-tab {
                   background-color:#eaf0f0 !important;  
                   border: 1px solid whitesmoke !important;
                   box-shadow: 1px 1px 1px 1px gray !important;
                }
            </style>"""

if __name__ != '__main__':
    light_style
    dark_style