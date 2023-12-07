#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 18 14:45:21 2023

@author: Dr Simon Lam
"""

from dash import html
Div = html.Div

def initiate() -> Div:
    legal = Div(children=[
        Div(style={'margin': 'auto', 'width': '980px', },
            children=[
                html.Hr(),
                html.H1(
                    'Privacy policy'
                    ),
                html.P(
                    'Last updated: 18-May-2023'
                    ),
                html.P(
                    'This privacy policy describes Our policies and procedures on the collection, use, and disclosure of Your information when You use the Service and tells You about Your privacy rights and how the law protects You.'
                    ),
                html.P(
                    'We use Your Personal Data to provide and improve the Service. By using the Service, You agree to the collection and use of information in accordance with this privacy policy.'
                    ),
                html.H1(
                    'Interpretation and definitions'
                    ),
                html.H2(
                    'Interpretation'
                    ),
                html.P(
                    'The words of which the initial letter is capitalised have meaning defined under the following conditions. The following definitions shall have the same meaning regardless of whether they appear in singular or in plural.'
                    ),
                html.H2(
                    'Definitions'
                    ),
                html.P(
                    'For the purposes of this privacy policy:'
                    ),
                html.Ul(
                    children=[
                        html.Li(
                            children=[
                                html.Strong(
                                    'Account'
                                    ),
                                ' means a unique account created for You to access our Service or parts of our Service'
                                ]
                            ),
                        html.Li(
                            children=[
                                html.Strong(
                                    'Affiliate'
                                    ),
                                ' means an entity that controls, is controlled by, or is under common control with a party, where "control" means ownership of 50% or more of the shares, equity interest, or other securities entitled to vote for election of directors or other managing authority.'
                                ]
                            ),
                        html.Li(
                            children=[
                                html.Strong(
                                    'Laboratory'
                                    ),
                                ' (referred to as either "the Laboratory", "We", "Us", or "Our" in this Agreement) refers to the Steve Jackson Laboratory, Cancer Research UK Cambridge Institute, Robinson Way, Cambridge, CB2 0RE, United Kingdom.'
                                ]
                            ),
                       html.Li(
                            children=[
                                html.Strong(
                                    'Device'
                                    ),
                                ' means any device that can access the Service such as a computer, a mobile phone, or a digital tablet.'
                                ]
                            ),
                        html.Li(
                            children=[
                                html.Strong(
                                    'Personal Data'
                                    ),
                                ' is any information that relates to an identified or identifiable individual.'
                                ]
                            ),
                        html.Li(
                            children=[
                                html.Strong(
                                    'Service'
                                    ),
                                ' refers to the Website.'
                                ]
                            ),
                        html.Li(
                            children=[
                                html.Strong(
                                    'Service Provider'
                                    ),
                                ' means any natural or legal person who processes the data on behalf of the Laboratory. It refers to third-party companies or individuals employed by the Laboratory to facilitate the Service, to provide the Service on behalf of the Laboratory, to perform services related to the Service, or to assist the Laboratory in analysing how the Service is used.'
                                ]
                            ),
                        html.Li(
                            children=[
                                html.Strong(
                                    'Usage Data'
                                    ),
                                ' refers to data collected automatically, either generated by the use of the Service or from the Service infrastructure itself (for example, the duration of a page visit).'
                                ]
                            ),
                        html.Li(
                            children=[
                                html.Strong(
                                    'Website'
                                    ),
                                ' refers to Deoxyribonucleic Acid Damage Response Clustered Regularly Interspaced Short Palindromic Repeats Screen Data Explorer, also known as DDRcs, accessible from ',
                                html.A(
                                    'https://sjlab.cruk.cam.ac.uk/ddrcs/',
                                    href='https://sjlab.cruk.cam.ac.uk/ddrcs/',
                                    rel='external nofollow noopener',
                                    target='_blank'),
                                '.'
                                ]
                            ),
                        html.Li(
                            children=[
                                html.Strong(
                                    'You'
                                    ),
                                ' means the individual accessing or using the Service; or the company, or other legal entity on behalf of which such individual is accessing or using the Service, as applicable.'
                                ]
                            ),
                        ],
                    ),
                html.H1(
                    'Collecting and using Your Personal Data'
                    ),
                html.H2(
                    'Types of data collected'
                    ),
                html.Ul(
                    children=[
                        html.Li(
                            'Usage Data'
                            )
                        ]
                    ),
                html.H3(
                    'Usage Data'
                    ),
                html.P(
                    'Usage Data is collected automatically when using the Service.'
                    ),
                html.P(
                    'Usage Data may include information such as Your Device\'s Internet Protocol (IP) address, browser type, browser version, the pages of our Service that You visit, the time and date of Your visit, the time spent on those pages, unique device identifiers, and other diagnostic data'
                    ),
                html.P(
                    'When You access the Service by or through a mobile device, We may collect certain information automatically, including, but not limited to, the type of mobile device You use, Your mobile device unique ID, the IP address of Your mobile device, Your mobile operating system, the type of mobile Internet browser You use, unique device identifiers, and other diagnostic data.'
                    ),
                html.P(
                    'We may also collect information that Your browser sends whenever You visit our Service or when You access the Service by or through a mobile device.'
                    ),
                html.H3(
                    'Tracking technologies and cookies'
                    ),
                html.P(
                    'We currently do not use persistent cookies, session cookies, or web beacons on our Service.'
                    ),
                html.H2(
                    'Use of Your Personal Data'
                    ),
                html.P(
                    'The Laboratory may use Personal Data for the following purposes:'
                    ),
                html.Ul(
                    children=[
                        html.Li(
                            children=[
                                html.Strong(
                                    'To provide and maintain Our Service'
                                    ),
                                ', including to monitor the usage of Our Service'
                                ]
                            ),
                        html.Li(
                            children=[
                                html.Strong(
                                    'To manage Your requests'),
                                ': To attend and manage Your requests to Us.']
                            ),
                        html.Li(
                            children=[
                                html.Strong(
                                    'For other purposes'),
                                ': We may use Your information for other purposes, such as data analysis, identifying usage trends, determining the effectiveness of our promotional campaigns, and to evaluate and improve our Service, products, services, marketing, and your experience.'
                                ]
                            )
                        ]
                    ),
                html.P(
                    'We may share Your personal information in the following situations:'
                    ),
                html.Ul(
                    children=[
                        html.Li(
                            children=[
                                html.Strong(
                                    'With Service Providers'
                                    ),
                                ': We may share Your personal information with Service Providers to monitor and analyse the use of our Service.'
                                ]
                            ),
                        html.Li(
                            children=[
                                html.Strong(
                                    'With Affiliates'
                                    ),
                                ': We may share Your information with Our Affiliates, in which case we will require those Affiliates to honour this privacy policy.'
                                ]
                            ),
                        html.Li(
                            children=[
                                html.Strong(
                                    'With Your consent'
                                    ),
                                ': We may disclose Your personal information for any other purpose with Your consent.'
                                ]
                            ),
                        ]
                    ),
                html.H2(
                    'Retention of Your Personal Data'
                    ),
                html.P(
                    'The Laboratory will retain Your Personal Data only for as long as is necessary for the purposes set out in this privacy policy. We will retain and use Your Personal Data to the extent necessary to comply with our legal obligations (for example, if we are required to retain your data to comply with applicable laws), resolve disputes, and enforce our legal agreements and policies.'
                    ),
                html.P(
                    'The Laboratory will also retain Usage Data for internal analysis purposes. Usage Data is generally retained for a shorter period of time, except when this data is used to strengthen the security or to improve the functionality of Our Service, or We are legally obligated to retain this data for longer time periods.'
                    ),
                html.H2(
                    'Transfer of Your Personal Data'
                    ),
                html.P(
                    'Your information, including Personal Data, is processed at the Laboratory\'s operating offices and in any other places where the parties involved in the processing are located. It means that this information may be transferred to - and maintained on - computers located outside of Your state, province, or country, or other governmental jurisdiction where the data protection laws may differ from those from Your jurisdiction.'
                    ),
                html.P(
                    'Your consent to this privacy policy followed by Your submission of such information represents Your agreement to that transfer.'
                    ),
                html.P(
                    'The Laboratory will take all steps reasonably necessary to ensure that Your data is treated securely and in accordance with this privacy policy and no transfer of Your Personal Data will take place to an organisation or country unless there are adequate controls in place including the security of Your data and other personal information.'
                    ),
                html.H2(
                    'Delete Your Personal Data'
                    ),
                html.P(
                    'You have the right to delete or request that We assist in deleting the Personal Data that We have collected about You. Please note, however, that We may need to retain certain information when we have a legal obligation or lawful basis to do so.'
                    ),
                html.H2(
                    'Disclosure of Your Personal Data'
                    ),
                html.H3(
                    'Business transactions'
                    ),
                html.P(
                    'If the Laboratory is involved in a merger or acquisition, Your Personal Data may be transferred. We will provide notice before Your Personal Data is transferred and becomes subject to a different privacy policy.'
                    ),
                html.H3(
                    'Law enforcement'
                    ),
                html.P(
                    'Under certain circumstances, the Laboratory may be required to disclose Your Personal Data if required to do so by law or in response to valid requests by public authorities (e.g., a court or a government agency.'
                    ),
                html.H3(
                    'Other legal requirements'
                    ),
                html.P(
                    'The Laboratory may disclose Your Personal Data in the good faith that such action is necessary to:'
                    ),
                html.Ul(
                    children=[
                        html.Li(
                            'comply with a legal obligation'
                            ),
                        html.Li(
                            'protect and defend the rights or property of the Laboratory'
                            ),
                        html.Li(
                            'prevent or investigate possible wrongdoing in connection with the Service'
                            ),
                        html.Li(
                            'protect the personal safety of users of the Service or the public'
                            ),
                        html.Li(
                            'protect against legal liability.')
                        ]
                    ),
                html.H2(
                    'Security of Your Personal Data'
                    ),
                html.P(
                    'The security of Your Personal Data is important to Us, but remember that no method of transmission over the Internet, or method of electronic storage is 100% secure. While We strive to use community accepted means to protect Your Personal Data, We cannot guarantee its absolute security.'
                    ),
                html.H1(
                    'Children\'s privacy'
                    ),
                html.P(
                    'Our Service does not address anyone under the age of 13. We do not knowingly collect personally identifiable information from anyone under the age of 13. If You are a parent or guardian and You are aware that Your child has provided Us with Personal Data, please contact Us. If We become aware that We have collected Personal Data from anyone under the age of 13 without verification of parental consent, We take steps to remove that information from Our servers.'
                    ),
                html.P(
                    'If We need to rely on consent as a legal basis for processing Your information and Your country requires consent from a parent, We may require Your parent\'s consent before We collect and use that information.'
                    ),
                html.H1(
                    'Links to other websites'
                    ),
                html.P(
                    'Our Service may contain links to other websites that are not operated by Us. If You click on a third-party link, You will be directed to that third party\'s site. We strongly advise You to review the privacy policy of every site You visit.'
                    ),
                html.P(
                    'We have no control over and assume no responsibility for the content, privacy policies, or practices of any third party sites or services.'
                    ),
                html.H1(
                    'Changes to this privacy policy'
                    ),
                html.P(
                    'We may update Our privacy policy at any time without warning. We will notify You of any changes by posting the new privacy policy on this page.'
                    ),
                html.P(
                    'We will update the "Last updated" date at the top of this privacy policy.'
                    ),
                html.H1(
                    'Contact Us'
                    ),
                html.P(
                    'If you have any questions about this privacy policy, You can contact us:'
                    ),
                html.Ul(
                    children=[
                        html.Li(
                            children=[
                                'By email: ',
                                html.A(
                                    'sl681@cam.ac.uk',
                                    href='mailto:sl681@cam.ac.uk'
                                    )
                                ]
                            )
                        ]
                    )
                ]
            
            ),
        html.Hr(),
        Div(style={'margin': 'auto', 'width': '980px', },
            children=[
                html.H1(
                    'Terms and conditions'
                    ),
                html.P(
                    'Last updated: 18-May-2023'
                    ),
                html.P(
                    'Please read these terms and conditions carefully before using Our Service.'
                    ),
                html.H1(
                    'Interpretation and definitions'
                    ),
                html.H2(
                    'Interpretation'
                    ),
                html.P(
                    'The words of which the initial letter is capitalised have meanings defined under the following conditions. the following definitions shall have the same meaning regardless of whether they appear in singular or in plural.'
                    ),
                html.H2(
                    'Definitions'
                    ),
                html.P(
                    'For the purposes of these terms and conditions:'
                    ),
                html.Ul(
                    children=[
                        html.Li(
                            children=[
                                html.Strong(
                                    'Affiliate'
                                    ),
                                ' means an entity that controls, is controlled by, or is under common control with a party, where "control" means ownership of 50% or more of the shares, equity interest, or other securities entitled to vote for election of directors or other managing authority.'
                                ]
                            ),
                        html.Li(
                            children=[
                                html.Strong(
                                    'Country'
                                    ),
                                ' refers to the United Kingdom of Great Britain and Northern Ireland']
                            ),
                        html.Li(
                            children=[
                                html.Strong(
                                    'Laboratory'
                                    ),
                                ' (referred to as either "the Laboratory", "We", "Us", or "Our" in this Agreement) refers to the Steve Jackson Laboratory, Cancer Research UK Cambridge Institute, Robinson Way, Cambridge, CB2 0RE, United Kingdom.'
                                ]
                            ),
                        html.Li(
                            children=[
                                html.Strong(
                                    'Device'
                                    ),
                                ' means any device that can access the Service such as a computer, a mobile phone, or a digital tablet.'
                                ]
                            ),
                        html.Li(
                            children=[
                                html.Strong(
                                    'Service'
                                    ),
                                ' refers to the Website.'
                                ]
                            ),
                        html.Li(
                            children=[
                                html.Strong(
                                    'Terms and Conditions'
                                    ),
                                ' (also referred as "Terms") mean these Terms and Conditions that form the entire agreement between You and the Laboratory regarding the use of the Service.'
                                ]
                            ),
                        html.Li(
                            children=[
                                html.Strong(
                                    'Third-party social media service'
                                    ),
                                ' means any services or content (including data, information, products, or services) provided by a third party that may be displayed, included, or made available by the Service.'
                                ]
                            ),
                        html.Li(
                            children=[
                                html.Strong(
                                    'Website'
                                    ),
                                ' refers to Deoxyribonucleic Acid Damage Clustered Regularly Interspaced Short Palindromic Repeats Screen Data Explorer (also known as DDRcs), accessible from ',
                                html.A(
                                    'https://sjlab.cruk.cam.ac.uk/ddrcs/',
                                    href='https://sjlab.cruk.cam.ac.uk/ddrcs',
                                    rel='external nofollow nooperner',
                                    target='_blank'
                                    ),
                                '.'
                                ]
                            ),
                        html.Li(
                            children=[
                                html.Strong(
                                    'You'
                                    ),
                                ' means the individual accessing or using the Service, or the company, or other legal entity on behalf of which such individual is accessing or using the Service, as applicable.'
                                ]
                            )
                        ]
                    ),
                html.H1(
                    'Acknowledgement'
                    ),
                html.P(
                    'These are the Terms and Conditions governing the use of this Service and the agreement that operates between You and the Laboratory. These Terms and Conditions set out the rights and obligations of all users regarding the use of the Service.'
                    ),
                html.P(
                    'Your access to and use of the Service is conditioned on Your acceptance of and compliance with these Terms and Conditions. These Terms and Conditions apply to all visitors, users, and others who access or use the Service.'
                    ),
                html.P(
                    'By accessing or using the Service, You agree to be bound by these Terms and Conditions. If You disagree with any part of these Terms and Conditions, then You may not access the Service.'
                    ),
                html.P(
                    'Your access to and use of the Service is also conditioned on Your acceptance of and compliance with the privacy policy of the Laboratory. Our privacy policy describes Our policies and procedures on the collection, use, and disclosure of Your personal information when You use the Service or the Website and tells You about Your privacy rights and how the law protects You. Please read Our privacy policy carefully before using Our Service.'
                    ),
                html.H1(
                    'Links to other websites'
                    ),
                html.P(
                    'Our Service may contain links to third-party websites or services that are not owned or controlled by the Laboratory.'
                    ),
                html.P(
                    'The Laboratory has no control over, and assumes no responsibility for, the content, privacy policies, or practices of any third-party websites or services. You further acknowledge and agree that the Laboratory shall not be responsible or liable, directly or indirectly, for any damage or loss caused or alleged to be caused by or in connection with the use of or reliance on any such content, goods, or services available on or through any such websites or services.'
                    ),
                html.P(
                    'We strongly advise You to read the terms and conditions and privacy policies of any third-party websites or services that You visit.'
                    ),
                html.H1(
                    'Termination'
                    ),
                html.P(
                    'We may terminate or suspend Your access immediately, without prior notice or liability, for any reason whatsoever, including without limitation if You breach these Terms and Conditions.'
                    ),
                html.P(
                    'Upon termination, Your right to use the Service will cease immediately.'
                    ),
                html.H1(
                    'Limitation of liability'
                    ),
                html.P(
                    'Notwithstanding any damages that You might incur, the entire liability of the Laboratory and any of its suppliers under any provision of these Terms and Your exclusive remedy for all of the foregoing shall be limited to the amount actually paid by You through the Service or 0 USD if You have not purchased anything through the Service.'
                    ),
                html.P(
                    'To the maximum extent permitted by applicable law, in no event shall the Company or its suppliers be liable for any special, incidental, indirect, or consequential damages whatsoever (including, but not limited to, damages for loss of profits, loss of data or other information, for business interruption, for personal injury, loss of privacy arising out of in any way related to the use of or inability to use the Service, third-party software and/or third-party hardware used with the Service, or otherwise in connection with any provision of these Terms), even if the Laboratory or any supplier has been advised of the possibility of such damages and even if the remedy fails of its essential purpose.'
                    ),
                html.P(
                    'Some jurisdictions do not allow the exclusion of implied warranties or limitation of liability for incidental or consequential damages, which means that some of the above limitations may not apply. In these jurisdictions, each party\'s liability will be limited to the greatest amount permitted by law.'
                    ),
                html.H1(
                    '"AS IS" and "AS AVAILABLE" disclaimer'
                    ),
                html.P(
                    'The Services is provided to You "AS IS" and "AS AVAILABLE" and with all faults and defects without warranty of any kind. To the maximum extent permitted under applicable law, the Laboratory, on its own behalf and on behalf of its Affiliates and its and their respective licensors and service providers, expressly disclaims all warranties, whether express, implies, statutory, or otherwise, with respect to the Service, including all implied warranties of merchantability, fitness for a particular purpose, title, and non-infringement, and warranties that may arise out of course of dealing, course of performance, usage, or trade practice. Without limitation to the foregoing, the Laboratory provides no warranty or undertaking, and makes no representation of any kind that the Service will meet Your requirements, achieve any intended results, be compatible or work with any other software, applications, systems, or services, operate without interruption, meet any performance or reliability standards, or be error free or that any errors or defects can or will be corrected.'
                    ),
                html.P(
                    'Without limiting the foregoing, neither the Laboratory nor any of the Laboratory\'s providers make any representation or warranty of any kind, express or implied: (i) as to the operation or availability of the Service, or the information, content, and materials or products included thereon; (ii) that the Service will be uninterrupted or error-free; (iii) as to the accuracy, reliability, or currency of any information or content provided through the Service; or (iv) that the Service, its servers, the content, or e-mails sent from or on behalf of the Laboratory are free of viruses, scripts, trojan horses, worms, malware, timebombs, or other harmful components.'
                    ),
                html.P(
                    'Some jurisdictions do not allow the exclusion of certain types of warranties or limitations on applicable statutory rights of a consumer, so some or all of the above exclusions and limitations may not apply to You. But in such a case the exclusions and limitations set force in this section shall be applied to the greatest extent enforceable under applicable law.'
                    ),
                html.H1(
                    'Governing law'
                    ),
                html.P(
                    'The laws of the Country, excluding its conflicts of law rules, shall govern these Terms and Your use of the Service. Your use of the Service may also be subject to other local, state, national, or international laws.'
                    ),
                html.H1(
                    'Disputes resolution'
                    ),
                html.P(
                    'If You have any concern or dispute about the Service, You agree to first try to resolve the dispute informally by contacting the Laboratory.'
                    ),
                html.H1(
                    'For European Union (EU) users'
                    ),
                html.P(
                    'If You are a European Union consumer, you will benefit from any mandatory provisions of the law of the country in which you are resident.'
                    ),
                html.H1(
                    'United States legal compliance'
                    ),
                html.P(
                    'You represent and warrant that (i) You are not located in a country that is subject to the United States government embargo, or that has been designated by the United States government as a "terrorist-supporting" country, and (ii) You are not listed on any United States government list of prohibited or restricted parties.'
                    ),
                html.H1(
                    'Severability and waiver'
                    ),
                html.H2(
                    'Severability'
                    ),
                html.P(
                    'If any provision of these Terms is held to be unenforceable or invalid, such provision will be changed and interpreted to accomplish the objectives of such provision to the greatest extent possible under applicable law and the remaining provisions will continue in full force and effect.'
                    ),
                html.H2(
                    'Waiver'
                    ),
                html.P(
                    'Except as provided herein, the failure to exercise a right or require performance of an obligation under these Terms shall not affect a party\'s ability to exercise such right or require such performance at any time thereafter nor shall the waiver of a breach constitute a waiver of any subsequence breach.'
                    ),
                html.H1(
                    'Translation interpretation'
                    ),
                html.P(
                    'These Terms and Conditions may have been translated if We have made them available to You on our Service. You agree that the original English text shall prevail in the case of a dispute.'
                    ),
                html.H1(
                    'Changes to these Terms and Conditions'
                    ),
                html.P(
                    'We reserve the right, at Our sole discretion, to modify or replace these Terms at any time without warning.'
                    ),
                html.P(
                    'By continuing to access or use Our Service after any revision become effective, You agree to be bound by the revised Terms. If You do not agree to the new Terms, in whole or in part, then You must stop using the Website and Service.'
                    ),
                html.H1(
                    'Contact Us'
                    ),
                html.P(
                    'If you have any questions about these Terms and Conditions, You can contact Us:'
                    ),
                html.Ul(
                    children=[
                        html.Li(
                            children=[
                                'By email: ',
                                html.A('sl681@cam.ac.uk',
                                       href='mailto:sl681@cam.ac.uk'),
                                '.'
                                ]
                            )
                        ]
                    )
                ]
            ),
        html.Hr(),
        Div(style={'margin': 'auto', 'width': '980px', },
            children=[
                html.H1(
                    'Disclaimer'
                    ),
                html.P(
                    'Last updated: 18-May-2023'
                    ),
                html.H1(
                    'Interpretation and definitions'
                    ),
                html.H2(
                    'Interpretation'
                    ),
                html.P(
                    'The words of which the initial letter is capitalised have meanings defined under the following conditions. The following definitions shall have the same meaning regardless of whether they appear in singular or in plural.'
                    ),
                html.H2(
                    'Definitions'
                    ),
                html.P(
                    'For the purposes of this disclaimer:'
                    ),
                html.Ul(
                    children=[
                        html.Li(
                            children=[
                                html.Strong('Laboratory'),
                                ' (referred to as either "the Laboratory", "We", "Us", or "Our" in this disclaimer) refers to the Steve Jackson Laboratory, Cancer Research UK Cambridge Institute, Robinson Way, Cambridge, CB2 0RE, United Kingdom.'
                                ]
                            ),
                        html.Li(
                            children=[
                                html.Strong('Service'),
                                ' refers to the Website.'
                                ]
                            ),
                        html.Li(
                            children=[
                                html.Strong('You'),
                                ' means the individual accessing the Service, or the company, or other legal entity on behalf of which such individual is accessing or using the Service, as applicable.'
                                ]
                            ),
                        html.Li(
                            children=[
                                html.Strong('Website'),
                                ' refers to Deoxyribonucleic Acid Damage Clustered Regularly Interspaced Short Palindromic Repeats Screen Data Explorer (also known as DDRcs), accessible from ',
                                html.A('https://sjlab.cruk.cam.ac.uk/ddrcs/',
                                       href='https://sjlab.cruk.cam.ac.uk/ddrcs/',
                                       rel='external nofollow noopener',
                                       target='_blank')
                                ]
                            )
                        ]
                    ),
                html.H1(
                    'Disclaimer'
                    ),
                html.P(
                    'The information contained on the Service is for general information purposes only. For research use only. Not for use in diagnostic procedures.'
                    ),
                html.P(
                    'The Laboratory assumes no responsibility for errors or omissions in the contents of the Service.'
                    ),
                html.P(
                    'In no event shall the Laboratory be liable for any special, direct, indirect, consequential, or incidental damages or any damages whatsoever, whether in an action of contract, negligence or other tort, arising out of or in connection with the use of the Service or the contents of the Service. The Laboratory reserves the right to make additions, deletions, or modifications to the contents on the Service at any time without prior notice.'
                    ),
                html.P(
                    'The Laboratory does not warrant that the Service is free of viruses or other harmful components.'
                    ),
                html.H1(
                    'External links disclaimer'
                    ),
                html.P(
                    'The Service may contain links to external websites that are not provided or maintained by or in any way affiliated with the Laboratory.'
                    ),
                html.P(
                    'Please note that the Laboratory does not guarantee the accuracy, relevance, timeliness, or completeness of any information on these external websites.'
                    ),
                html.H1(
                    'Errors and omissions disclaimer'
                    ),
                html.P(
                    'The information given by the Service is for general guidance on matters of interest only. Even if the Laboratory takes every precaution to ensure that the content of the Service is both current and accurate, errors can occur. Plus, given the changing nature of laws, rules and regulations, there may be delays, omissions or inaccuracies in the information contained on the Service.'
                    ),
                html.P(
                    'The Laboratory is not responsible for any errors or omissions, or for the results obtained from the use of this information.'
                    ),
                html.H1(
                    'Fair use disclaimer'
                    ),
                html.P(
                    'The Laboratory may use copyrighted material which has not always been specifically authorized by the copyright owner. The Laboratory is making such material available for criticism, comment, news reporting, teaching, scholarship, or research.'
                    ),
                html.P(
                    'The Laboratory believes this constitutes a "fair use" of any such copyrighted material as provided for in section 107 of the United States Copyright law.'
                    ),
                html.H1(
                    'Views expressed disclaimer'
                    ),
                html.P(
                    'The Service may contain views and opinions which are those of the authors and do not necessarily reflect the official policy or position of any other author, agency, organization, employer or company, including the Laboratory.'
                    ),
                html.H1(
                    'No responsibility disclaimer'
                    ),
                html.P(
                    'The information on the Service is provided with the understanding that the Laboratory is not herein engaged in rendering legal, accounting, tax, or other professional advice and services. As such, it should not be used as a substitute for consultation with professional accounting, tax, legal or other competent advisers.'
                    ),
                html.P(
                    'In no event shall the Laboratory or its suppliers be liable for any special, incidental, indirect, or consequential damages whatsoever arising out of or in connection with your access or use or inability to access or use the Service.'
                    ),
                html.H1(
                    '"Use at Your own risk" disclaimer'
                    ),
                html.P(
                    'All information in the Service is provided "as is", with no guarantee of completeness, accuracy, timeliness or of the results obtained from the use of this information, and without warranty of any kind, express or implied, including, but not limited to warranties of performance, merchantability and fitness for a particular purpose.'
                    ),
                html.P(
                    'The Laboratory will not be liable to You or anyone else for any decision made or action taken in reliance on the information given by the Service or for any consequential, special or similar damages, even if advised of the possibility of such damages.'
                    ),
                html.H1(
                    'Contact Us'
                    ),
                html.P(
                    'If You have any questions about this Disclaimer, You can contact Us:'
                    ),
                html.Ul(
                    children=[
                        html.Li(
                            children=[
                                'By email: ',
                                html.A('sl681@cam.ac.uk',
                                       href="mailto:sl681@cam.ac.uk"),
                                '.']
                            )
                        ]
                    )
                ]
            ),
        html.Hr()
        ]
    )
    return legal
