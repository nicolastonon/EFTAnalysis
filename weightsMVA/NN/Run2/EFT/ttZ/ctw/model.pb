
A
MYINPUTPlaceholder*
dtype0*
shape:���������
S
Feature_normalization/CastCastMYINPUT*

DstT0*

SrcT0*
Truncate( 
�
Feature_normalization/sub/yConst*
dtype0*U
valueLBJ"@   �J`@   @�L�?   ���h@     �]@    $���   ��e@   ��(Q@   @Lȑ?
b
Feature_normalization/subSubFeature_normalization/CastFeature_normalization/sub/y*
T0
�
Feature_normalization/truediv/yConst*
dtype0*U
valueLBJ"@   .�t@   @J��?   ��w@    'Os@   ���@   @/-w@   ���a@  ���3 @
m
Feature_normalization/truedivRealDivFeature_normalization/subFeature_normalization/truediv/y*
T0
�
dense/kernelConst*
dtype0*�
value�B�2"���������5�Ⱦ���n���=Gٽ^C�?�mf>��	�Ͳ�>A��=D!>�U>���e?�d�>]Ս��ߩ=nm?_�I?���>J��=C�����?	f���>�9�>��u�S�E?=l\��+?�`?�:���Z)?u?�Z��\�>B7&��JJ>�ǚ?7l}>7m�>_���>z�H<�7�t����Rؾ.r?�����>�?�Y=��=7�?�ِ?�?���H(?��r���?;E�>in�>���K�>�܄��5��ヾח1=JT㾐JE=��J�������<���X�K�3?�Ef�����%h��H�?�Q?�L��}4?c�"?�'X���������O<�R?�/��nK?/8�>�h;�:��>�l�pw�>�5 >��e?_y� X?7�s?�	?|{A��b!��2�=�\F=��ʽD&=��+��^!=E�=�����,��u�Z)%=��?���ܾ���=���>�	�������lV>6�>t§�i���оE����^�>�.1>��,�?F;#�&<�,��>|�>��T?ژҾ�`�>��T�/�a>Ej�>Yݠ��ھ�'�>�R����G2>�1f)?J�>s� >7?��p�[��}G�>6>�}?��+��D�?�ק�V��*���d
)���c�ǜs?hRe>Ԯ�f����OP�;�پ��3���*>/?���?M��wl�=�:`?êQ?,�Ծ��p������#�=o �>�]�?T�Z>�j^>u��H?���=���>E���4�k���FϽؙ߼lZ�����6�����\?��:�g��=���f�;��Cվ��~��B">n�>��6=�r[��,���?ZA?t�s?�{>�y��.ɽ$V�=�.y>�"B?��>�l�_�Ӽy�I�-=�����=�(�=�����Dj?
5>��?we�=7M��{�Ҿ��k?�>��6*`���*?���>^��>�b?�'>?1��m׍?�c?>*�`?��Z���,���>q���MW��x��Jd�wQ�oԘ�`�>1jV�^u>���ad?��8?b���_�	��O,����>��>�W�6�Ҿ��m��΄��B�mk?��'�̏s>2��>����C*>���=��>*��>��>fF?��}�c�Y?ߟ�>��?&�P�@�=�h�{��?1�w�3�2��鑾|S�>�Ƥ>B�1>���>`��<*����.��5S?p�	�6����.��b>��?�I�=g��>,�?��?���S��>L��>�4?�;�=2V���f��[�B���><}m��~�,�$C��ݗ��@ W?����>W�z> iB?���=�&���ҾN �ǅ2?%}��"-�_�=���^?8­>�E���.�>P��k�<���?F��>� ���W���F?����a��>�N �E�=�]㾔4�>�H�b ���?�>��<�TZ=�o7�3#�=��ž�A�=o��W�7?v�(?�<�!�>�8��t�>��D>Ze��l��>L2�3��>'��L>�{�<�6־t�[����>G�3�L�w=s=ʕ?�t���׽>��>��A�m;�=ӧ"?��-���>R�?5v'?h��
�

dense/biasConst*
dtype0*�
value�B�2"����<2�=��;���=�3�=:���5�<F��47�����i���j�pM~<��D>zE�=bS�������½�s>�V,<�1$�.x��8��<"=�����4 �C�,<TA��A�ֽ;C�n"4=.q$��=�����4Ǝ��1�йھ�_C����� ;��i����ݘ:E�>9Pu=����o`�=��<K���w͆�
Y

dense/CastCastFeature_normalization/truediv*

DstT0*

SrcT0*
Truncate( 
>
dense/MatMul/ReadVariableOpIdentitydense/kernel*
T0
n
dense/MatMulMatMul
dense/Castdense/MatMul/ReadVariableOp*
T0*
transpose_a( *
transpose_b( 
=
dense/BiasAdd/ReadVariableOpIdentity
dense/bias*
T0
d
dense/BiasAddBiasAdddense/MatMuldense/BiasAdd/ReadVariableOp*
T0*
data_formatNHWC
*

dense/ReluReludense/BiasAdd*
T0
�
batch_normalization/gammaConst*
dtype0*�
value�B�2"��Cz?��u?��^?�{?[?�+�?�3�?�s?S6�?6	�?t�x?� �?ߣk?�in?֎z?��f?�?���?$Cr?a�V?щe?��?�HY?au�?`�o?a�O?��Y?_�b?ٗ?���?�$|?,ۅ?C΅?ַx?��?�g�?o�|?9:�?F�?��j?_Y�?`��?�L�?�%|?V�m?N�?O�S?��h?#�?��m?
�
batch_normalization/betaConst*
dtype0*�
value�B�2"�'�7>�,�LE�>����ő(�N�6>��<�򒼍��M��>��f�dZ��ښ>��6����>-�e������=x�;=5�6Y���ς>R����z.��n��3r7=,Ʈ�w54�3SǾN��>߽ز���d�=�X.��If�?���ҫ!=u=�>N�>H����+��	�=�*��m5�>�+�G돽�[>�́=t[��
�
batch_normalization/moving_meanConst*
dtype0*�
value�B�2"����>�nj>y�\>��u>%j >�ń>�;�>�>&�:^�5>���=H�>�OS>j��>]!>:k�>�[3=�S>Y)�>�Ď>$�=_��=L�J>f��=�v>>%!=b>&Ҋ>M�D>�=���>(�'<��=U?�>e�>#�w>ۖ=�,>o��=�y?>�g�=Xe|>:�>>�>��>`�<x�>�q>ed�> �>
�
#batch_normalization/moving_varianceConst*
dtype0*�
value�B�2"�c�>ced>�N�>�Y�>�O>r��>�_�>�Z>��:d��>Z!>?�F>j�f>D�)?�0>�?���=h�>��>��>�p�=���=�t�>zS�=R�>Lt�< ؏>/��>U�k>D��=v�@?I�R<zNp>�ۺ>��?�1�>���<%5_>s��=7Fz>�>��>u�@>
P�=-Z�>k<�_?���>���>�X2?
f
,batch_normalization/batchnorm/ReadVariableOpIdentity#batch_normalization/moving_variance*
T0
P
#batch_normalization/batchnorm/add/yConst*
dtype0*
valueB
 *o�:
�
!batch_normalization/batchnorm/addAddV2,batch_normalization/batchnorm/ReadVariableOp#batch_normalization/batchnorm/add/y*
T0
X
#batch_normalization/batchnorm/RsqrtRsqrt!batch_normalization/batchnorm/add*
T0
`
0batch_normalization/batchnorm/mul/ReadVariableOpIdentitybatch_normalization/gamma*
T0
�
!batch_normalization/batchnorm/mulMul#batch_normalization/batchnorm/Rsqrt0batch_normalization/batchnorm/mul/ReadVariableOp*
T0
b
#batch_normalization/batchnorm/mul_1Mul
dense/Relu!batch_normalization/batchnorm/mul*
T0
d
.batch_normalization/batchnorm/ReadVariableOp_1Identitybatch_normalization/moving_mean*
T0
�
#batch_normalization/batchnorm/mul_2Mul.batch_normalization/batchnorm/ReadVariableOp_1!batch_normalization/batchnorm/mul*
T0
]
.batch_normalization/batchnorm/ReadVariableOp_2Identitybatch_normalization/beta*
T0
�
!batch_normalization/batchnorm/subSub.batch_normalization/batchnorm/ReadVariableOp_2#batch_normalization/batchnorm/mul_2*
T0
}
#batch_normalization/batchnorm/add_1AddV2#batch_normalization/batchnorm/mul_1!batch_normalization/batchnorm/sub*
T0
�N
dense_1/kernelConst*
dtype0*�N
value�NB�N22"�Nt�d�=Y�[>#�v=���1.X��h�<�����!�=�,!=mM4�e�ջ���>:⸼��
��P�=��1�Y��rX=�(;>��;��EX>��|�d��������.���ܽ=v#>���B�����־?�%>�QH<��>c��>�� =~u�P `>E�Q�\�V�a=�<�d<�����ýʧ��\�}�>�c�����1Y��Ӓ>���=�A>��/�xݣ;d(
�\*A����o�\>���=<�>>=��c��&���\�q�d��T?9��=3t�>7��>����èS�3��=���]�c�Y_~�^g7="��;���<M��=�	���a>E��n�������*>0ZI�:]e������'>��=>X3��\B>��+<>��܅������x��L����<�L>�Ծ*j>~�ս>�9=UG>�)-�:̎>z�>I�/z��L��C�ƾ
��KPL��R���ch��	>Z��>����>���=R�H�׎�<��̾���7d�>����=y#�F1��U��=����-�<S�>[9
=�u��f�2��=�Q1�xL��G�>��锽�_=V����j��Q9���);��s�>��*>�>��о���Z������O־y��=l>�d!=�cj=�#>�o�<A�w>���>�<\��"������=�'�>`iV>��'�ܒ<�'�{lv��']=!�v��� >���n#�>��Ὥw#�0��>��\>w8O=��l���M��wu�(�����B�/��`*��u�:ՆA��3�<�[��A]���#	�ڵ>����Ե����>7��H�)=L�<>L� �O����?x���@���]�%��Ft��bu,>K� >5��>��^�;�4>�5�>>�>��x#��_*=��>u���<�F�Z���!��L>��bq@>�޽���A^f�a` >9�s=!�-�T����<��>��>�Ee>��>�c>/>�>��
����7^���=IIH>�� �8�j�T���Ϝ	<�u���2/>��=�U����<r@����#���V�&)ջ�]Ż��><Y�>�/��+��7��V� �T�^��k�*�>��%<�K�"�2��+�6::=*ej����<��>�<��˼�j�=>���=X���g�=fwټ+�$��	c>l�`��	>�x�!ؽ�$>�s����Q��t?����A�>	�'�>N���7�>{��`EZ>�]��N}���;��gk۽#���!������#c=H5X���>"᡼���>�~>��?~{'�!�}=��2>���`}�qi}>�R�=BE����k���_�����þiP>���^�e=�p>�>�f>{�=
�\������!������4���%�;�>Ľs>�'�>?�����>J��6��>��o��>��=2�����佋E;↾釖>y���B��=-�6>�-t���0>p��=����g���+���=D悽�~_��!�=&z�>а#>尪>�?9=�f�>��>��$>q��
؏>V�y>��O=~��KK�>S�`�mZ���u>]f]>��=òI>0 >��U=��>��i=�Ѕ>�1�U�ɽ�~C�?��>@����j?�/>��<}�>E��>ч�>�;&?�V�> 9	> ��=���=�I�=/�=q�>/H�>Ӽ�����B�=6H�>Hn�������>n@x>�:�>]�)>��� �Žm�=��3>l��>[��&�E>}l�>����(|��"!}���>���><��=/w�>���>�@��nv=5��>-	U��^>�@�>t�l�~.V��>�z�=�A����� �>hs��M:�=�q��(	7��~��M[��觢����=>B��~����N>kOV�k*3>`�,�n��\�<I��>�0]�i���=ꉽ�Y:�Z�>yt����\x���,���>���=62�<ڸ;�^���Cʾz�����*<�>�m)��۾���>����=3#�>x�<A%F�Ѵ��w��>�?�>��=�b(>�++>�����=�?-�c�>�)>��1���-�Nr>z�</���><?����>�˥;$4d>��@>�%����?�(=�>�B��IN>�!�x�>�8��U�=��i�>�c=�L׽���>_u7��g�=�!j>�b7��&��XOӾF�1�>e-�.n�ќx>�S��s���_>�J���ͯ���&=߫s=qy�>��,�.�ս"������>H��>�(>wBL=�P ���[>ԁ?�W{��c<�@�����=�����>N(�=z��<ѣ?vm�=���<�_W>܅�V�D��>-ࢾ�
�>������>�eo>���=�)��6r�=C�ڽ�2��M�XG�=�1B�^��>��[=�cZ>dF�P��=���i�=�>t<�=uȬ<H���� ��S�T��!�?��~��>����\e��ݻ��c�>ؓ��j =�3Ľ�S`=�nT��F���c=?AO>TK����z�'�m�=�2�><�=w�*�Ώ:>�"5��~Z=�^�=Ӗ%��?�=��5�n����Q=瓣>L����]><p�@w��oa?>:]j������ު�$��s����KS��E>���=@��=�1Q>����qx�s���0��>7�=FFZ>$3��Uƾ����k������$��=���<O�Q���>���=?�q>*]�2�[�Mh#=m�ļ;.s=�7���>����f>��>9���8b���EL)�a�?���=ڤ[=1٠��dy>�=�D=�\��,���]��� � >�+�<�-��6���ĕ�:�>���>k�p�>O�M�0����+��<��#�iI��M���쨾� ����������d�ѽ�=�=��,���>�&��Σ���Ҿ��߽�s�>!��=��C<�e�� >�e�>��I������x��R>CW`=�킾R��ɔ+������'��T$>�|�C�c�����Q�S�%���Ⱦ ����k����<�V	��M>�ǽ����`��}Q��H'>(F�=#�������O�<�+����=Y�>�Q\����=\��F�>�c�G�軙��>�]F>�达Ѕ¼1�>oE�=2�8�H#D��8��^�н ��=�n�;2V�=P�8>��R>��N>ؘN�Ñ��޺�>��M�; >=�hJ>h	Z����>9�>�>)��>oVG>1�P�Q��=<\���6�>��>���>�<4>�������j���b��AH����.=�S�>
H>���>���7��>����!��*>� �2��B����t@��A�=G.>Jq�>:
R>���Ͻn�c�P<5w0�6#���[�>͜ټ��>@5(�s��]7�=��h���>��Z�Q��>��>u$!��.>�eս��߾�>�z
>/�X>���K���z�����-���#��6`��S2�sZ�zT!��,<=�d��">��>~�~=ֽ���M�l�:�<a;3>��d��A�>�=�>����k����˾@9��n%z���>s��g�>C;��7Y�<W�:�,QO>��S>T6[�ϰ<�<lQ��8a�Z� >�{�=9��=�W��f��f�#�<���ؽ9�Ѿ�Ew�Mq���dŽ�+=�����g�/s��.���d�>U�a������E=O�1>�dy�JB>&3>���mK�=��o>��.=OJ���-*�|��Ǩ�=��y=�����\6=�`�=L�>M�սl�=!,R=ݺO���o��:;�h5�>Aߡ>��%>�ld>/m<>�^�>���
��1=t1��A���c�=�{���>���="
��wJ�֍�4B>>|��hk�=҇�>u9w�4߃>�"Y���=P�>���9�=�l�����<�����I���ꓽzk��k���Y0m��6�>rk�=�oǺ1W�	S���\D��X9��,˽���\����Y>�5<nя>��>J+=>���=�=��,���k�>1N7>pyJ>��"��k�<����&�=�ao>鵏���<�U`=Pm>8G>c�U�T���5n־���>�~@>|�;ʄX>��>DU>ջ7>P��=pB�>�͟>R�?�]5�=U���v�轍�{>�b>y�����>���;?w�'�R=w~���;��>7]Ժ��H�'d=s�-������p��T=�Mȋ;�E�=ͥŽA��=�|1>/�;>~d=��W�=%;=�>4(�}��� پz̾>1վY�+�D+�<l�}=@{��$C��́��]=R@�C���RV�=)�u��Ǿdt�PLu���S�bp=6ι<'�V��Ǡ� V��/�꽁��=�;׽րy<Ѥ�>f�k��;�>��о	�=rX��Kn��$�=b?=%�1�C��>��GL��5g��պ=�?>J<����=��3���=n�I=��?>��<')�<%�=�r8> �:� ��>3%�<�5�����<�F>�fB>"�>�;#=? ��~׻��Y>в�>m�>��=�P��&�C=e&>7h�Ǝ>}�+Xb>Wٽvt��'8C����=���'�ѽ��Q������vg�͸ �KE���|�=*c>���=J�%�{W��n�=u� �f��f��r�>q	��^��=.����A=�n��SY�<�-������䝾�>�eL��ǂo��f.>�m���ʃ��Y>��>e�-�s����|�>"b	������3��	�6��=P����>ޯ��M�{>� �>ަ��b���;ó>���=p+���B�>���:�N�=`�?=���k�=��~>`��>����YJ�>ɮż�S���^�m�>�薾�:V��qi���Q<�6�>�=���>c�����#�>��=�b���M��pVy>�8�>4��<ѐ�s����9;��>���=RĎ�ae)>��s������>���>�w�>X�m��c?��/n=��Y>��-<=�>��-<۾�'��,4��e3�<�O���N>���a>r="
��S��>�����>,��<<>�|t�տ�=j��=��X<��2���[=�憾��y>0�žB;�>�H�:��>Ň��Q�㽂�>8�j!> 5c��=�dּ�U.�o.�=���� IU>�.-�Z�>nz=�޽� N�a���Ya= �	��;>�d����;�w��h6�>uŽQ=�ΫO=�j<%�I�t׶�/8�=�^->���k���F�<Q{a�MB}>����7�B�	�h�$�~�:��	K>d����a#�[!�=_bz<%߽5�!��<+<�14���=+a>
ъ=�4>��>N>]?"=<t=�	�=�>,�Ծ������ =�l=ҏҽAn>D��>�����@��i.=�z߽������>~�7�r�h>-����y>ԿG>oͅ��G�����=iˇ=���>�,�=��>}�U��;�=��b�O����d>�J��P����s>q�=���N̾���>�E��_��>�+>�5>�<�6�=8}>�w�=s�Q>�j��f���z�BB>��=��_>��y=���>��D�����}�1>S^��K��!�4=M�_��.�����=�a��&�j�=.�	�F>�̊<���
�>˷���.���Q�W���y�T�d>���>7S4��6��>9'�Ia����Ὦ�=�˻Ȇ�bh���ܳ�t
=X>Z>W}>�½&k�X�=��<>������p{ý�>d��Ũ�ϫ��'�4�>s[>G�B=�|=�e�>J侧آ<`f��O�Y��>xȄ� x�����zڽ�ƻ>MS�CU��R�<݈,��"��LXV�J��=� `����=q���
>���>'{�>�n?>�Rʾe	�/-��������>!	�>��ؽ�0�=&s�����3o>�5�>�>ʹ������=ac>��'�
)�=xa?>�u�=��>4-���������ûV��>M�i<g4}�@��= ?�]�p�G^>]�L>}��>�?d�ʾ�!R>�r>��_X[>W����ꂾ�d����>���>�S[=ݩ��k"y��E��O�>Ύ`�ϼ��\�ǽCz��RJL��cO�I�< ߝ>��9�hh}��7�=��>��>5����q�<�>Z���~��;�i����T�<,$�=稲�֗ؽ�>k3���\+>)��>��=Vj�<L�[>3��<�8\=�Ĉ������2�>�_�>3w�>ex�s��vG�<��	>_na�K�\[3=)w���~ν�� �I�+>��> �>����<�u@>�M�=��)>��>P�<>��w��KV>b��=b\@>�5�@�<w�>��=�=���<��>�9A=�3?��=c?ǽ�{�=,6�=Uk��tڅ>�.㼒(����mN�=�)��D�>b=����>�X����=�x>�W��1�=��v>��a��Դ�hM	��I=�ӽ�cG�#���9����ߣ�=j:~=� >�j�>���Pmɽ�Æ>�<������ߓd>�O���&�=W��,w>%p�G�+��Q�=��&>��=�Y>8H�:�雾L���T�=�8�=�2>�x��<DJ>��>�+�<YA/���>�o_>���>:'����=/�Q>*V$�'Au>��=���=��P�~mϽ6�x=mV���+�숾���i�=P��=������h�p�����E����e��c�>Th>^����j>s�����Ծ}h�=J����
=|X���
3�j�6<<'�;�꾵)ýk�=���>	@>h1=:�\=x?���E"������˞>L�<@:p>��b=�eW�U��.I��c�����^rI��(=�5=�^�<�2(>�.�H\�>�6|����>�\���>Ep�>��>R0���˽!�=_9��ҵ=���>��$������>��>ͧ���/����5>5�v>�7�-O���f�c(���������>��d�]\�=U����)�����>�#�>�Fξ�~>��C�a*c>l)�=q^��`���>��Ľ�\>�P+�$�q>�Y��mo���Q�Pi>y�=屆>f����>w� >l�>K΍>�C���˽��$>��>�@>D�=���*ߚ�m�=>�{>�}�W���pm���<ɺ$������%�>�� >���>�3!<q��;�:>m3�v���i��{��#�>Q��;a��潰��LA>!��=�(�=>;�}�=P(>�+��.T>hb�=]@��_w��S�<%y)=�>;��=A\I>#z�=���>���*����=V>=�>9ُ=3�=h�B>��=��T�B>.a�>���{��>k��&J���j�� �<�ۍ=e�=+�'=?�C��>'�>��>�.e>�~��}M>ҹ>��>��>R�/���!�z�����*=�@�Ӯ>W ���^�>E>>=f�>p� ��	�v�>�	�>���ar>�{�<r���R����ꔽ����BC�=� ��*;ų���l�=�V�vg�=oH%����&؍� ��;�ʹ>�?�>�P>jZ��׳�%`�5���Ԫ>��:a��֝B�	$о<�>���E�ܾSue;��y>k0�=	�w>��=d��X�����Q�3���-��=�*������@����_&�������>���=����E���=�i���)�=��W��LW����#�W~�����	��N1I�A=Ľ�i��Vގ>���=B��<�D!���]�ȋ¾P�*��m�P��S>��Ӿ��O=x��>�{�>k���a߭>�`���)��ﯾ�Y���=>l�ýP` �;������p����b���c^�Ԥ)>I�#��І��K>�w�>%d�<�#�7�S��>��#��=U�Q���Sq>���UƐ=?2����F�8�޾���h >��>U�>� ��o������=䊲��~�����N��� �6��oSѾ��`���>�%����=� 4�xz>2`��_��c��E�xڨ=�2'>>�>�ᓼ�86���=�G
>Ȩ8���>�{��<̾7�>@vR�/>k�i,�=7���jO=71�=�7=�i�q=M�>ܴľ���>�-=f����>e����W>7;>~[T;��>�N�=�H[�4�>�^X>�P%>�?�>�;�=��!>�r���6�=��>,3@���>��ν�N����.�V.�>�_3=9?�����ެ1<47R=��]>���,س>�f�=��	�L��_�g��XQ;��>��k; &�
T�>�Kn>����w�>�'">�yZ>��|>��Y�ܨe=DI�����%ʩ>n�G=���>�]�>T<z>������>zY�=�M=ɗ	>�ʦ�if�>>��@��=�g>��
>l�>G?U�W�ܪ��x��<���=S+
=G��>�ҫ�ij>\���c#>FUj=�Ǧ>����.�x>J�λΒ�<.�<=��{>/�>�N<��>�M��˱��cɾ��ʾ� �>�� >( <O�>���:���1	��t��Z�/>�G��*�ܾ�&����+�/�ܾ���}�P��>��4=v~0����j���x�����=��&��n6�h��>ͿǾ@U���þGlN=�����>�T��u���᰾����h�J����>��-�k�������7����>'6>\@�>9ힾ�d>Xd��㶛�x˰�'=B��龃9a=�Ѥ=Aɐ=�t�>��%>(�'>AK=d�W>�\=�}u>=�@�$w�=��^������z��G�Մм�h�>z��<վ��)ʻu<"�Ƥ�<y�O>T��5-����H�H�=/_�����>�ۼ�B��cr�5~5��������Y��L >�Z3=՟ȼ�YԾ��e>�搽���������P>�K�<5��>�T��R�U�����������k2��_Fb�K�?>���=��:�b7Q��&,>�!��G���ok�M1�����=]��}�:��y>�%�>��վ��F��$˾=è���g��=��=W�X>O&�ǚ~���s>%�y�ܐ>7z��T=��>2L��3f>��L;U��=yѽLw���;�>\�V=�:L<p3�>���>@X�܅�>�N>�*=ieD>�{�>�)=��=���P��=�}	�7v���l>en�=�l6>"	8�bN1>UF����>�>�e���h�>�儾Wʢ=%�t>d<���H�չ�{>�X��2�=#Y�>x��>�b>�)����$ο���o���=s��=��ʾZ����$:�����$���'�>�<���=b�>	 W��S�>��t=����YR������9X>����.���jD��8=<�>������]�he�rj����>���>���(oؽX`h�H��>�a��=�<�>-{b>ߕ=Q"X>���=v��=�������>Hm��n<�8>k�g=-ي��uE=��>�3�<�����"i=;�`��*��w>������>���=9�
>褝��V=��u���a>N{�<欯�;�'>�iG��D�+Tо��#����=%�ѽ���8�D>�\���Cd�����=�����=?R��ލ�>����B�=�L=S�B��4ݽ�j޾��K�#>��ܒ�>ۗs�	4����> s�=M{���8=鎕��ᅽ�N$>����=�!��+�(>i��=cXӽnؾp?,�ڸS<��A>�%���P��^�VMb=5��������=ԇ7>�l>�n���sW=h>KuZ�GD)>��>. �z|�>g�>�Y>��k���=��6���=�IS=�,�=Q�*>�g�������j �=��>WҺ��9�t�ֽ��x�b��=&w>��}�����
j�����k�>X~ >��=u=��罴2>>Sھ>W��(���v�j�N#��2�~���>T��>羘>�>3��ʀB�W�:C>>D}�=�O�>g=�q�=�>�U��� �>��Y=s��<x�>��ؾ���=۳�.�^�l��=Oּ>����	��I>!�;5�ü,>>T��;��Z>b��>�-G>/��>\|p���0>���=������W;
�
dense_1/biasConst*
dtype0*�
value�B�2"�v=4���Ֆ�|N�=b�uG!��у=���=!N�0�=�m5�������P��!�����u�E�;qi	�����O���=꾽_�<��Ⱦ��^�O����<��1��p��Z<V0�<�gR�1�I�8�t��H�=�4���o��A���ｲB�������MBi��'���q����;�y��S׼�=�F�
B
dense_1/MatMul/ReadVariableOpIdentitydense_1/kernel*
T0
�
dense_1/MatMulMatMul#batch_normalization/batchnorm/add_1dense_1/MatMul/ReadVariableOp*
T0*
transpose_a( *
transpose_b( 
A
dense_1/BiasAdd/ReadVariableOpIdentitydense_1/bias*
T0
j
dense_1/BiasAddBiasAdddense_1/MatMuldense_1/BiasAdd/ReadVariableOp*
T0*
data_formatNHWC
.
dense_1/ReluReludense_1/BiasAdd*
T0
�
batch_normalization_1/gammaConst*
dtype0*�
value�B�2"�uto?},n?��q?�`l?*�~?{G?�X�?��?��d?��]?2�u?��?~��?,��?�ڈ?(e?�ye?oa�?��?Vom?k[�?�=�?�C�?z�i?�BJ?:��?α�?K�?��L?�6\?�qj?�9�?�h?��u?��?GE�?qh?��k?-�\?��f?uwn?'��?%eC?��?2�?�8{?H�{?=�?��I?���?
�
batch_normalization_1/betaConst*
dtype0*�
value�B�2"�x��>����tv���=y"P>�w������ �M����=����@�D=j�>�P�=^r>�#.>;�u�s���d8�x�>*��9�<v�q>M��>fޤ�H����I�&��>;�J�U��=`�q�dլ=�>�5���)����;��=��=3<Y�q���������=w�>�/�=��>E�-�m����c�iRн��?��D"�
�
!batch_normalization_1/moving_meanConst*
dtype0*�
value�B�2"��cS?��?Ez�=�/T?�o>2
�>	�=?�9?c��>Z�>�>�p�>.[�>���>���>��?�3?�>x>��>S+?E�)?5Q�>O��=8�?9��>�~�>P�|>���>��?��=?�*b?�Fm>���>�v>�Ns?�g>�xy>?�X?2?�(�>�ץ>�zu>�3?��5>?�|>�"�?K(�>��>�q�?��?
�
%batch_normalization_1/moving_varianceConst*
dtype0*�
value�B�2"���M@)&�?�,u?�5�@0[0@�_C@d=H@vX*@��7@9�?s�?X� @$$;@p�O@��@�d@Q@�%@tb+@���?�O@I�@w7H?���?��?�;@©+@��:@���@I9`@��@q��?`s�?�8�?̠p@X8t?~��?�@��?��p@�?�K�@?@Vy�?��I@�@�:@��?}�X@��@
j
.batch_normalization_1/batchnorm/ReadVariableOpIdentity%batch_normalization_1/moving_variance*
T0
R
%batch_normalization_1/batchnorm/add/yConst*
dtype0*
valueB
 *o�:
�
#batch_normalization_1/batchnorm/addAddV2.batch_normalization_1/batchnorm/ReadVariableOp%batch_normalization_1/batchnorm/add/y*
T0
\
%batch_normalization_1/batchnorm/RsqrtRsqrt#batch_normalization_1/batchnorm/add*
T0
d
2batch_normalization_1/batchnorm/mul/ReadVariableOpIdentitybatch_normalization_1/gamma*
T0
�
#batch_normalization_1/batchnorm/mulMul%batch_normalization_1/batchnorm/Rsqrt2batch_normalization_1/batchnorm/mul/ReadVariableOp*
T0
h
%batch_normalization_1/batchnorm/mul_1Muldense_1/Relu#batch_normalization_1/batchnorm/mul*
T0
h
0batch_normalization_1/batchnorm/ReadVariableOp_1Identity!batch_normalization_1/moving_mean*
T0
�
%batch_normalization_1/batchnorm/mul_2Mul0batch_normalization_1/batchnorm/ReadVariableOp_1#batch_normalization_1/batchnorm/mul*
T0
a
0batch_normalization_1/batchnorm/ReadVariableOp_2Identitybatch_normalization_1/beta*
T0
�
#batch_normalization_1/batchnorm/subSub0batch_normalization_1/batchnorm/ReadVariableOp_2%batch_normalization_1/batchnorm/mul_2*
T0
�
%batch_normalization_1/batchnorm/add_1AddV2%batch_normalization_1/batchnorm/mul_1#batch_normalization_1/batchnorm/sub*
T0
�N
dense_2/kernelConst*
dtype0*�N
value�NB�N22"�N!t>L�o���=�#=,;�=e�d>M>q{����f<���<���� !>���<�ɿ>�q�;u4>�\�B�T�={rZ=h>)B>~>2��n|x=>��>+�
>xm����]�	)=��K��$��a�!=�b{!=� Լ^Ib���>�X�>YT�=��>#�d>�>�&�>4��=1���0ˡ�:)b>
�f�[�ʽ���������2�<X{2�}��=�7�]4>����>��x��<���>�s�>��|�����e�S���%>-��r�ڽ󍌽�&�=zh���>hZ}>n��>@U��O->~Ҿ��3=�v���B�q=�ϑ� ���S	��,��=:'P=V�>�1��J�#>��1���s�������<�7�=T�y��=����E�W��(�>(������:FѼ��~= ��>0���l���b�>��S>ꟳ����.
�>=ؐ�<�C =ۦ�<��Y;���8F>�}^<x���FG?�������H�"e����	?�rĽ܂=��%>��&���N>4+�����=u�W�Tb�=!U>�=l=�=	��v��'�=\�>2� ��p�=%
?�Y�=ԨG���ݾ��6���&=c�v}�����f�>��
>Z%�t�>�����ٓ�E_�>s4->)��>|�h>��-��Ӱ��=�1�������p>Ѫ)>��=���>+м0���vI���w_5=��=�!�:�A�>�>P�\>l��յ�=������y��Z<�U���=��Ҽ~�u>�6�ᳱ�;񽉀�=��=i��>Պ�>@8�<�H�S#�>?��t/Ƽ)�>ʐ�>�(>��>�'c>e�
���־)��>>�|=Fg�Ey=ID�9�>�B>�_������D�>��N�ÂQ��=��>X�ؽi��=.&x>��)>�2>C��q��>�]�>U�}��)�����5V>~�I>	�L���>7��<�{�����=8u���Vn��_�>�ܶ=��9�4���<@^���i������U:($>�ʴ��\_�H��/='�R�˻ұ>h���������@n�����=G������>:��=��N���q>>1�=+�}=>˜�:��>��6=��Z��5=>H�Kŭ�i�>S3D�td<<�%3>YJh�`�>+溾��<@ ��w�>4����_�>�%\����>�7=Z�>5q��<"X����bE,����;L��L���`�����>�� ��U�>�%>�F�>s��=D(��+>�Ċ��:Q�ih�>���h<������-��=�mM>����ܠ=3h����>zs�=�tͽ��x>�魾�Yҽʂ��_����C3�ٿ]>b��>��j>��g>���>��>�4H��4��Tp�-��=�h>��>������9����>�ӷ����=h��K���_��=7��=�� =~g⽑51<�[�=���>��W�_�Z>Z�7�'���g�U�2>�h�>p	5�O��}v��I��p�=��<oǾ~��<{�>�����"=pM���_�> �F<T��*�`���댖�T��>ҫ�=6X�>�,W>�B�լ���pw�1�>BX�>��>�a��@�^�"�k>��c<��7��==~k���>Zz=o�*=��>���n6���>�;��܏%��/˽��=(��>���mu�����;�N6�+�o=��d��uo��E><��>[�->^�<K��>��>��g���F=!��=�u�O�����=��z>��w>/c�>�3�=`��=���>˷�>x��=�sM>�s|>b��<�Dl>﮾>�|N>�c�>�b���|�(b;=�p��{�ʼ}���A��>_�=����x2S�j�⯽3��>g�9�gԖ>l>ԟ��f>��ʾ3��>�����+��������G�>R�V>7ؽ/i0=g� �%��<�V>Cn>k�A<��5� <�<:XK��v>�=�S����W=��[�Mt�P�R=�+߽1�>1s��c%�{[�=Y�,>�7W>�|�=��>�A�>oZ�Dk�<��v��>`̙�]�����>�$/<���E�^>d M�wGK>��+";�/3��qw�3�%>��>O�0>0ֽm~�>V)9�5y:�@?�>� m�z�>�W��{�=�g��Q���&Y<�����=�)�jz
>���>�`=P#>���;����E�B�)>G��<��<>� >M3��/2�n?��J�>=������=�@�>:A���`����">[q���d�2�\>.�k�_��<�ȁ>P��;"c�y.ѽ����1�:>pⱽE�>�h�<����=W����=�_9>0�C=�3+>���q�<f�>�uν��=e��^<�>O�ͽ�=A�>R�½�/���S�>v��>+\�=�=�ީ���<��y����=��+=-�}���=�0�<�<�t?�{ƾE�=�*���c�>we�><�>�Q=i\��L�= �������>�F���J5���q<�=v�4>��+���rg��M�4�pjW>�8���y�=Q���j_�=H��>�K=n�����=�
˼n"��V	O>��<��=���>�'8>ɕ=P|Z��,��+�(>ľ��2�ge��r���8��>��+���R��4>%�?��o>�6>�g�&i�>�Q�>葽�[��=6(>Н㽝��=<�>#�>�E��cr�J�,��H��W����
>a�H��J�<�6޽i��;�J�TBҽ��;f1�~��>�P�=�l�>��ĽzK��%΄>s�u=Y������=\�+�����a�X��ꪾ��<,���꘾��:��Հ>n,�>AZ =�B�>[�o���?�=�q��<���(o>
�>�(�����>P@���W=Ď>�>8>������,x�������O�X�<��Q>[�2��9��?�q��Qf<��o<a �=H�>Y���yP<�/?x�v>i}q��L�=�)0�1��>�n>B�>E�����|>���>���>(�_��=��Q>٧�=��>+m���.�<9���� ���>�\����<0�J>*	^>`ۏ�D+B�Rд�z���sȌ�N�����P��>�j+>Jޏ>�仾��>]���CQ��=?= �p���߽z�w��>�>%��>[�N>U�@>��������˾���>o��>�7��̫T<���\}��Uc>���<2Ҡ�;�i<pw���1>ńY<=�ɾ��#�L���VTν���>ܵ�l,=���<�@ ���=�ͽ�̟��:j�;��4��Ͽ^>�az>W��>�U�Ø־�D�>��>~���I��"Ό=ȵ��\��=b��=q�=8?=f�>懽�m����>4����#�Ϧ���>��D=���`؊=#4�f0�=�񁽩�K=�%�o����Q�>�?����>ڀ���5���Ao��A>˖��wr�a#�>aa �����=a��\�� e>ʝ�<�"��-�a����>O�\>'݇��3�<oe;���>�f�(�F�>�	K����=�_�xz��dG@����=�����&>bГ=�uh>�Co�4x>���vk�>ǧU��g=�E���=��8>��K>��񽫳o=~y�>�f\��(����1�
S��+7�z��=����'���c=�)�=E�=pe�>>���Z�?�)��n��=�\:�U�u<�DȽ`F�>��&��4<$�>�'e����<,Ý���<��H>Y>�>�G�< Ѥ�ꥇ=��R�?�>�I�>i�c�5X6>�1�>$ʽtck�D�'>۲6>2�z��4�;m��7����>L��zY����9�����o��۷>N5R>j�<R"r����>�� �)�!V��v_>+|f���Ծ���~�g���ݾm�
=��b��ƾ=W�]>�B�!�̽O߾���~��6��5��>Sv�><�6��#$�C�<�L���~�s���&>,�x�\L�=���������Ƚ#>1> �*�V�����+U;�Ѱ?�%�>2z	>E�U=�� �=�$�%>��s�)����K>��>o��=q�����>?<?���:f�=�.U>�P�>~v����=w��u�<���;���<p6�>r�K��4>�á=y�>,F >)�u=lEp=��>��F�^��kB�(�l�2�W>�~V<�ǎ<0�B�qQo=�U%��\>!>b*����<�7��)>_%\����W���"-�DC>�r�>�F	�#K</�2�`d����=���0Ь>��b>g3�>W>!8�>�=�� =���q��J�$��E%c>8/�>XB�=�+�=����ٷ>ka7�b���`�v�b�$�H�;����Y����=6��>1q�=r7=5ú�hOe=V]���VW:X����?�=��=�J)���t�!uK��&�=��r>;3/>ޝ�v;�u��E��<c���>*�ѽ).>S����=Ҷ�=��O>[����9>m�B=�󛼘P �,=������˗��u�>��p=��`��o>�Y�>�?��ܓ3>�ܤ>�mj=3�佷�SxM>*��h����{�&��>B'���VW=�b�l/=?���{?ߋ�<mTH>���=H��>�d��$����=>����|�=,?Z=[}1>��>iR�=�2=�d�=��>�J�������H>��	�2<2����=���z���V�=�J�q�>�2>P�۽�MN>���ͻ�=��n��)���������;���ּ7n� �>j�l<��<��>�OU>�~߾B�D=0��=�5=��>�����Y��<���=���=���=#|нG�t�3X��NUz�=��>�]ʾls5=Y�$���>���Bw�>�,>�J���J ��W���Iͽ<�ľ
�ɾ�L��+Ӎ>�����`�2�^���S>Vin>��W>]�ɽ�\<�Ơ�]#����Ͼ���:�+>���B�>��=���T`ƽD>�='$�?�;:��j>[��=Ҳ�<[���矾7������.�'����=B��<�)g>Io�=Q�C>�>�}��V`��IO���=*g���O�/�M>GD�>Rw�>��=�-�L>0�Z��L[���оO�*��ܱ=GS�;�g>�x_>���=q�Z�?u �}��]$��c>�w���cf�{��=]�/�MR��t��|�>=0�>��/��=���=�(=r���>+o�=^�v�֌�=� ľ̫_>u��=�P�>c�&�V�r>71-=Ĉ�=���=�^�>��;����o�3> �쾓S�>@-+>��f>�}���:<��0�]��>�	>ω>���>�v�>�#�=Y�>Ӷ<�[$=w'���n����8>�.������Y�-�;vU���6ʾvM��?�7�IЇ�E[�=��>�e_>�Y���LG=�ֱ=��>��i#>��<z�>� >�E�>����W@�T>�<)Ϯ>�G���>$��=T�=�ξ>/>���=?m+�;�	�O��0��^ к퀅�j�{=�<
>lؾ�٩>�}�=���~�>�)>��Rss���aB�<���|�=\��Z�:�>�=�Ҵ�N�>;�����>�?j��$���)=��2�����:%�=�֏�����>?�.>}*>kY�=��>�y��n���=��{=�u>^�V>�3 ����	S�=¹=�\C���<�=�8�G��yF�> ~E�A����;�O�=t b�1ܲ�2�ľߺ�=vK���>f�(=JA�>�:�������/���=�э��7Ҿ'���wC�����=`��>�ʂ>@I8=Yr���ɽ-=�"�=ß9�k(�>G?��O ��>�v��>�)��m=�>Bv�>PY&�ѐ�='y��C�8����=�>j��L�Խc��=���e�>�q�>�~q>n�R>9�
>�=g��#�����<=C�<3	�!-�=���=�p�>m���$�;$u�>��<�1���:��Ln>� ��:��/f<�&���1�>��a��8�>u3>�M�=��#>����D2�� �/��=%�=AHF�SO�=]�J��@��>k�ɽ;?M��3���J���->��>d��;������]>q������ O���>�}W�7ˆ>���A�>�3���>�Lu�;��:�%�=S �<:C>�+�>��b>����1���+��q�>sw�=8�j>�ˏ=Q~��nн�����=JfC=D%��|D�W��:���U&�����>#Ġ=�Q~=Cc�=`*�>
�网�	�)�4>�.y>�C>K�ڽ5I���6Ӽ.��>�g�>�������D��v����=�<�>
s;��O>0�߾u��=�H�>$�}�Ȯ�=J/>�q>>]��>�d>Σ(>�19�^%�E]���s>�K�>������H>�<���~�N(>N{ۼ��>�?����??�¾��<S)Z�!`�>Hb��w����>=S$�>*J>�_��}X>)s�����c�F=�n�����пA>�H��w�v�\o/��нw�<.>oim�En%>�������,� ��=O;�=����Vk;�(9�"�����-��-�=��v>��?>��8����>K*>잹>A�O����~�<g%���J�6�>��Ⱦ�|�>�=eܻ=o��=ؼ�4����J2�:�a��>�"�^y�(v��5��=ަ��4U<�
���X�>>֟�Ҹ��xt?�B��ƣ>�^����>��=�b�c�l��Lo=�=>����0b�;sI��ٙ����>^?7�,7��i��@��T;>u�>�>��K=�V�<�������>e�C>��q��I�=�(�9�#���!��z�=)��=�h>[��>���>U���~F�<��o��m�>���k�>��X�֥{���<�>�9;�qѿ>,�evQ=��#> �G�YH�>�8_����˰��={���M>��[>�	M>&�H�����ɍM�֋4={�=tQ�=�^4>?��=�q�=>��>�Օ=F!��E�>!������>�=�x�,=OR�>k�Z��S#>���8�,��`����=Y���M��>�ٴ>߫�6 6���!���=B��Q��=��>M=H�o[�?�>�p=>S6�4;�>g>�9� 뤾�>ó���u���4'=_��>���g�Y���/e>ƙ\>\!k��R>���=z�=��<�1���V>�*Ծ�"a>�ّ�,��>�`e=@]>���b�>���]���Ɔ�g�l=T^0>-�}>\��N���p�����>8��=�̬>���"3��u�=��޽X��<Z��>%��_�~�,�=v�=Q
��6��������)�,�>���5�>#5�>�=���2 >�b>&��>��H�|��>�ڦ����=�Ѣ��.����u>
]�!d�+mB>��=&
���t���
<�|U��OX=ѽ!�8�T����>�S�=����x�<�$�<�+�=�I�=Y%)�X�\��/���%=�� ?d��=�*�=0W�>�^�=1�=�m�=�R>�������Ϋ�>?�>�[���5<���֑W�я�m�,>�C=�#�&���q<���>���<d�>X��������w���{�>
(>�_�c����4P=A"=�*�=W���Q�Kȭ�j��`8>�Iՙ������:�>U��So�0��=þZA�<-z>��@��?>2m�^b�>N���jQ�>�؂�����4=�<��K��&gվ0��x��.��W���x�� �=>�,G����4ʑ�N��;�ҏ��^>�BO>�"7��h>Vg���!>����=�v=5�X������D���r���v�U�J�r�M�yce�������{<}�����="Z7=#_�=m�߾�M��j\��r���d>�;F�y�?����=}�=���#��z��/x����=;譾Z�>�U���><(U>O=ž��=ʍ->��
�rK���l4>U=��
�δ@>�'=�eF���w�E�<ö�3���Ŋ��V�D>l�v<irS=�Y���RD>VT�<�q�=�Ž>�ѽ��,�E�Ks >&�����C���P��Ւ=���=E��f.���ҵ>^�~>Z19>C���Y��>l�i�nܿ��$@�8.ܽ��>ٟ=�w�������Bm>�WýC#Ƚ���c>���=t>-�\��a4>ؕ'�~�ɾ�5�=�!��cA>c�;�̼��;����=��|>yϞ��� �=�>o2>,>�>�>�N�	�Q>xO>��.�U!=�έ����+}B���ľS	>��q����t9�>(ؙ>dR�%@>���>�O>U�C=�!�q�ż\��E>IϾyz>��S>�|�=G�=�t�=�d>ɮ<񢤾65�=�j >@�%���/>���D�=�)�,�t>���>��J߸���>=OG� �����,Ȧ>����� ��fӽ���܊>��w����>ig>|��>����>��A>���<�-m>�fM>P<>�*�=/l�>r�ͽ��>/�<��+>��=��=kۈ<��`�Z��)x=�}$�2!�>ɿ���>��>�xY=N��A�=8)H>�ǈ>������>6� >,�<���ži��`i=����<юJ=y������N�p�>K�������,�R������>�]�`����<?0>�t2=�c�>Ⱦ����<k��>�s�>�[+>R��>�*�>�B�=^Й��s���zS>�N�=Sl�=��=(湽�/�=p��=`iC��&V�V?�i���>,|>��%<c����7�=^�R��;��*>�#>��=4j�'�>aE)=
aj�@��>��0=]��>c3%�,�=h*#>)ؼ��<*.��P"��?R�?Lq�<G�I>�M��L><8q>֚>��9<�>&���k�þ�/Ҿ��=R<���p9>�^������"	�	����(�t�W�t�%>y�<!b>�н7AL>E����V`>�����A��=̽=��\x�=�����>	.>��ھc��=G�Խnn�e+>�U�¼?�4�Q4�>*e�����<Zȫ�#�=�����>��ڽ�<}�1D��a���F>�i�=���=q�V�$=�����7�b�K=��l�K�=�,>���=�[
�t2"��~`>�>�>�>^�k�Օ�<�p'>����2ω>��\>��=cG�z��=���>G�>pf���(�=�`�=ʖ�=��>!տ<<��=�x�>�n���>��߯��=�m�=ke�=�=J����Z>�U>+���,��>��S>�����f�>(��acC��rH��ج>�=�\~>�A����B�gߠ=�m�=@��b>m}�>0���X줾�g�s�\�$����=���v����>��e��>�c\�Eci�Q@0���2=mkM��6һ<N >�8���*>P�>^F��v�=��g>��>���>��O=9*�<�;����=�L�\#��Zٽ�>����=��M�D=[���>W&=S�'= ��>[A;��~���=P"�<���>(��</8%��T��.8��˴�R�=:z>x�>�`ս��h��;��`>2+4�x�=<�ᾅ��>>�v>Ы^=�&W>Ez�X�>W�9>���{�I<�Ȃ>/�0�|lý��Q�=�}�f0{�s5?_���Dxy>B'�J�������~K��N>c�><�!	_=�*����+>1�������s���]c}�U�F��x>��|>�ݟ>����'�<��F��)�7�.���޻��V����=߬2=���>X,���Ͻ�D>���8�P�#8׽G�q����x!>8�ѽ<�L��N߽0M+��x��M�>��=+C�yH���=R'>�>ↁ>O�1�r��X�=>���������g��)>Y�>IԼ�r�*���=K̾C��=�K7<��Q>S�=/d��{�k�%>Qd��.�=v��G�=馈>n�>:q=Z�ȽUE�S9�*�.�����+�x�x<����R>�Gt�}e�<#`���j	=�=<Q����g<?�6�>%��=�u�ӫK�����/�<�>M�l>��j�����
�
dense_2/biasConst*
dtype0*�
value�B�2"�� �=�GH>�~������z��v-���>���=.�<s�=�2ٻd=y(
>P�=�����/��<���=��=ԭ->-�P>�Q>3=C=���=}�]��U>
q��#l�="�>��<>ٹ�=\�j�ޅ����=���=�-�� I�=�'�=�8>>g->Ǿ�=�E>��=(>��b�6�Ϣ=H��=ּ>��(>
B
dense_2/MatMul/ReadVariableOpIdentitydense_2/kernel*
T0
�
dense_2/MatMulMatMul%batch_normalization_1/batchnorm/add_1dense_2/MatMul/ReadVariableOp*
T0*
transpose_a( *
transpose_b( 
A
dense_2/BiasAdd/ReadVariableOpIdentitydense_2/bias*
T0
j
dense_2/BiasAddBiasAdddense_2/MatMuldense_2/BiasAdd/ReadVariableOp*
T0*
data_formatNHWC
.
dense_2/ReluReludense_2/BiasAdd*
T0
�
batch_normalization_2/gammaConst*
dtype0*�
value�B�2"�G2>K1?b��>���>��(?L?x�>�*?\S*?�&?%�>��=�>V)�>?]��>���>���>~'�>�`�>��>���>�Î>�T�>�=�>�T"?,��>u�>ǐ�>D�>��>Me�>6�9?lŬ>�#?�G�>=ϟ>��>B�*?���>�>d�?���>�'?��m>/X�>s�>j2�>���>r�>
�
batch_normalization_2/betaConst*
dtype0*�
value�B�2"�8�����<B�u<��;�6��+��:=�T��W�;�Z��6�$<�!��f~<�e����32�<��Q�<XH<d�<{Y4��X���=7X�<
"��T�!<3Ԃ�y#O�drK�~��<:������=,{���?<1�;|�;����?`�ϛ�;��r�p�*�t�
��K�q�IA~�S�Q<kې�,$=�E<
�
!batch_normalization_2/moving_meanConst*
dtype0*�
value�B�2"��1�>���?��W>�)�>$�Z?n�n?��G?7�X?�v?8j`?�g�>�&?�Ta?ն?(�%?��>K�>�?)��>.*?KB?��*??p�5? ��>�C]?SOY>���>VG.?��?� ?I�>�Fs?�W0?W�7?d�>�Z)?Q6?A/o?V;U?2�?��E?�5?��t?x#�>W9�>iv2?�	?���?�{?
�
%batch_normalization_2/moving_varianceConst*
dtype0*�
value�B�2"�,?۩]@���>)B?��
A㮽@ yo?�ջ@�F@��#?a�x?t8�?��M?�8�?��A@2}�?V<?)�?{N�?dc?|d&?���?l?�?��D?��?��?�s?��>J�?��Q?]ڔ?��H?�c�@c�?�v?M�n?�K#@��`?���@�a?1)@;�?m5�?�/�@�R"?�k�?��?�B?��?���?
j
.batch_normalization_2/batchnorm/ReadVariableOpIdentity%batch_normalization_2/moving_variance*
T0
R
%batch_normalization_2/batchnorm/add/yConst*
dtype0*
valueB
 *o�:
�
#batch_normalization_2/batchnorm/addAddV2.batch_normalization_2/batchnorm/ReadVariableOp%batch_normalization_2/batchnorm/add/y*
T0
\
%batch_normalization_2/batchnorm/RsqrtRsqrt#batch_normalization_2/batchnorm/add*
T0
d
2batch_normalization_2/batchnorm/mul/ReadVariableOpIdentitybatch_normalization_2/gamma*
T0
�
#batch_normalization_2/batchnorm/mulMul%batch_normalization_2/batchnorm/Rsqrt2batch_normalization_2/batchnorm/mul/ReadVariableOp*
T0
h
%batch_normalization_2/batchnorm/mul_1Muldense_2/Relu#batch_normalization_2/batchnorm/mul*
T0
h
0batch_normalization_2/batchnorm/ReadVariableOp_1Identity!batch_normalization_2/moving_mean*
T0
�
%batch_normalization_2/batchnorm/mul_2Mul0batch_normalization_2/batchnorm/ReadVariableOp_1#batch_normalization_2/batchnorm/mul*
T0
a
0batch_normalization_2/batchnorm/ReadVariableOp_2Identitybatch_normalization_2/beta*
T0
�
#batch_normalization_2/batchnorm/subSub0batch_normalization_2/batchnorm/ReadVariableOp_2%batch_normalization_2/batchnorm/mul_2*
T0
�
%batch_normalization_2/batchnorm/add_1AddV2%batch_normalization_2/batchnorm/mul_1#batch_normalization_2/batchnorm/sub*
T0
�
MYOUTPUT/kernelConst*
dtype0*�
value�B�2"�T�)�Ғ;��5<9t�:���<ρ<i��;~'X<�ߺ��)<qɦ<��5<�<�&<*�
<lI���^Y<�m���;7!E<�sѷ���9�j�;]����(��}��qx�<�.;�J�+��;��\<�D�;Q�<L(�:y���a�$;�Y�;Ң8;mk�<,꛻#��:��C��<�	�<H<�<����D��;���m���
>
MYOUTPUT/biasConst*
dtype0*
valueB*j�b:
D
MYOUTPUT/MatMul/ReadVariableOpIdentityMYOUTPUT/kernel*
T0
�
MYOUTPUT/MatMulMatMul%batch_normalization_2/batchnorm/add_1MYOUTPUT/MatMul/ReadVariableOp*
T0*
transpose_a( *
transpose_b( 
C
MYOUTPUT/BiasAdd/ReadVariableOpIdentityMYOUTPUT/bias*
T0
m
MYOUTPUT/BiasAddBiasAddMYOUTPUT/MatMulMYOUTPUT/BiasAdd/ReadVariableOp*
T0*
data_formatNHWC
6
MYOUTPUT/SigmoidSigmoidMYOUTPUT/BiasAdd*
T0 " 